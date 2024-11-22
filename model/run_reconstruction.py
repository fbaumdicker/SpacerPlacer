import logging
import os
import pandas as pd
import copy
import pickle
import numpy as np
import json

from model import model_tools
from model.helpers import misc, import_data, stats
from model.bdm_likelihood_computation import symbolic_lh_computation
from model.experiments_tools import multiple_selection_fcts, construct_tree, filter_fct, \
    pooling_for_parameter_estimation, \
    tree_handling, name_inner_nodes_if_unnamed
from model.reconstruction_tree import ReconstructionTree
from model.data_classes.advanced_tree import AdvancedTree


TREE_GENERATION_PARAMETERS = {
    'simple_insertion_rate': 6522.48,
    'simple_deletion_rate': 137.20,
    'simple_alpha': 2.7278,
    'precise_insertion_rate': 0.6108,
    'precise_deletion_rate': 0.1830,
    'precise_alpha': 3.3377,
}

def run_reconstruction(rec_parameter_dict, dict_crispr_groups, save_path=None, plot_tree=True,
                       lh_fct=None,
                       logfile_path=None,
                       do_show=False,
                       sim_as_rec=False,
                       combine_non_unique_arrays=False,
                       tree_path=None,
                       tree_save_path=None,
                       use_provided_trees=True,
                       hide_unobserved_spacers=False,
                       selection_fct=None,
                       plot_order=True,
                       significance_level=0.05,
                       group_by=None,
                       finite_bl_at_root=True,
                       logger=None,
                       extend_branches=False,
                       tree_lh_fct=None,
                       tree_distance_function='likelihood',
                       tree_construction_method='upgma',
                       tree_gain_rate=None,
                       tree_loss_rate=None,
                       tree_alpha=None,
                       visualization_order='nb_leafs',
                       provided_aligned_arrays=None,
                       provided_dict_duplicated_spacers=None,
                       provided_numbering=None,
                       spacer_labels_num=True,
                       alpha_bias_correction=False,
                       rho_bias_correction=False,
                       core_genome_trees=False,
                       metadata=True,
                       alternative_parameter_estimation=False,
                       save_reconstructed_events=False,
                       dpi=90,
                       figsize_rec=(None, None, 'px'),
                       ):
    """

    :param figsize_rec:
    :param dpi:
    :param tree_construction_method:
    :param save_reconstructed_events:
    :param tree_lh_fct:
    :param seed:
    :param metadata:
    :param lh_fct:
    :param provided_aligned_arrays:
    :param provided_dict_duplicated_spacers:
    :param core_genome_trees:
    :param core_gene_trees:
    :param rho_bias_correction:
    :param alpha_bias_correction:
    :param use_provided_trees:
    :param tree_fmt:
    :param spacer_labels_num:
    :param provided_numbering: dictionary with spacer numbers and respective colors in tuples for each crispr_group_name
    fixes numbering and colors for reverse reconstructions (or in general)
    :param visualization_order: Just if custom order is supposed to be given for toy examples. Although I can implement
    other stuff in the future.
    :param tree_alpha:
    :param tree_loss_rate:
    :param tree_gain_rate:
    :param tree_distance_function:
    :param extend_branches:
    :param group_by:
    :param finite_bl_at_root:
    :param sim_as_rec:
    :param selection_fct:
    :param logger:
    :param significance_level:
    :param plot_order:
    :param hide_unobserved_spacers:
    :param tree_save_path:
    :param tree_path:
    :param combine_non_unique_arrays:
    :param rec_parameter_dict:
    :param dict_crispr_groups:
    :param save_path:
    :param plot_tree:
    :param logfile_path:
    :param do_show:
    :param alternative_parameter_estimation: based on reduced number of events (only considers spacers existing at root
    or only events that are in blocks with spacers existing at root).
    :return:
    """
    minimum_nb_of_arrays = 2
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    if logger is None:
        if logfile_path is None:
            logfile_path = os.path.join(save_path, '0_logger.log')
        logger = misc.create_logger('reconstruction', logging.INFO, outfile=logfile_path)
    logger.info('Running Reconstruction with parameter dictionary: %s', rec_parameter_dict)
    run_timer = misc.RunTimer()

    if lh_fct is not None:
        give_lh_fct = None if lh_fct == 'simple' else lh_fct
    else:
        give_lh_fct = rec_parameter_dict.get('give_lh_fct', None)

    if isinstance(give_lh_fct, str):
        logger.info(f'Loading lh function from {give_lh_fct}')

        # logger.info(f'Lambdifying lh fct...')
        # give_lh_fct = symbolic_lh_computation.load_lambdify_lh_fct(give_lh_fct, save_lambdified_lh_fct=True)

        # for later when I saved the lambdifyed fct one time
        give_lh_fct = symbolic_lh_computation.load_pickled_lambdified_lh_fct(give_lh_fct)
        logger.info(f'Loading time: {run_timer.time_from_last_checkpoint()}')

    ls_dict_protocol = []
    ls_boring_protocol = []
    ls_skipped_protocol = []

    if tree_path and use_provided_trees:
        logger.info(f'Loading trees from {tree_path} ...')
        file_ext = os.path.splitext(tree_path)[1]
        if file_ext in ['.pkl', '.pickle']:
            with open(tree_path, 'rb') as f:
                dict_trees = pickle.load(f)
        elif file_ext in ['.json', '.txt']:
            with open(tree_path, 'r') as f:
                dict_trees = json.load(f)
        else:
            raise logger.error(f'Unknown file extension of tree file: {file_ext} !')
    else:
        dict_trees = {}
    new_dict_trees = {}
    dict_data_for_lh_ratio = dict()
    new_dict_provided_numberings = dict()
    new_dict_provided_aligned_arrays = dict()
    new_dict_provided_duplicated_spacers = dict()
    for i, crispr_group in enumerate(dict_crispr_groups.values()):
        dict_bg_colors = dict()
        in_between_runs_time = run_timer.time_from_last_checkpoint()
        logger.info(f'Working on: {crispr_group.name}')
        if crispr_group.repeat:
            logger.info(f'Repeat: {crispr_group.repeat}')
        logger.info(f'Progress (group/total): {i + 1} / {len(dict_crispr_groups)}')

        skip, reason = multiple_selection_fcts(crispr_group,
                                               selection_criteria=selection_fct[0] if selection_fct else None,
                                               selection_parameters=selection_fct[1] if selection_fct else None)

        if skip:
            logger.warning('Skipped\n for reason: %s', reason)
            ls_skipped_protocol.append({'name': crispr_group.name, 'repeat': crispr_group.repeat,
                                        'reason': reason})
            continue

        ls_array_names = list(crispr_group.crispr_dict.keys())
        ls_arrays = [crispr_array.spacer_array for crispr_array in crispr_group.crispr_dict.values()]

        prev_len = len(ls_arrays)
        dict_combined_arrays = dict()

        if combine_non_unique_arrays:
            ls_arrays, ls_array_names, dict_combined_arrays = misc.remove_completely_same_arrays(ls_arrays,
                                                                                                 ls_array_names)
            nb_unique_spacer_arrays = len(ls_arrays)

            # logger.info(f'Arrays that were combined to one array {dict_combined_arrays}')
            logger.info(f'Note: For completely same arrays only one representative is used for reconstruction and '
                        f'visualization, due to provided option: --combine_non_unique_arrays. '
                        f'The classes of combined arrays are found in detailed results csv under '
                        f'"combined array names".')
            logger.info(f'Number of arrays before combination: {prev_len} /  after: {len(ls_arrays)}')
            # Be aware: Tree generation or alignment produces an error for single arrays (which are useless anyway)
            if len(ls_arrays) < minimum_nb_of_arrays:
                logger.warning('Skipped because arrays were combined. Too few arrays after combining non-uniques.\n')
                ls_skipped_protocol.append({'name': crispr_group.name, 'repeat': crispr_group.repeat,
                                            'reason': 'too few after combining uniques'})
                continue
        else:
            ls_unique_arrays, _, _ = misc.remove_completely_same_arrays(ls_arrays, ls_array_names)
            nb_unique_spacer_arrays = len(ls_unique_arrays)

        if dict_trees:
            logger.info(f'Loading tree from imported tree dictionary...')
            if crispr_group.name not in dict_trees:
                skip, reason = True, 'No tree provided for this group'
                logger.warning('Skipped\n for reason: %s', reason)
                ls_skipped_protocol.append({'name': crispr_group.name, 'repeat': crispr_group.repeat,
                                            'reason': reason})
                continue
            tree = import_data.load_single_tree_from_string(dict_trees[crispr_group.name] + '\n')
        else:
            if tree_lh_fct is None or tree_lh_fct == 'simple':
                modifier = 'simple_'
            else:
                modifier = 'precise_'

            if tree_gain_rate is None:
                tree_gain_rate = TREE_GENERATION_PARAMETERS[modifier + 'insertion_rate']
            if tree_loss_rate is None:
                tree_loss_rate = TREE_GENERATION_PARAMETERS[modifier + 'deletion_rate']
            if tree_alpha is None:
                tree_alpha = TREE_GENERATION_PARAMETERS[modifier + 'alpha']

            logger.info(f'Constructing tree with distance function: {tree_distance_function}, '
                        f'gain rate: {tree_gain_rate}, '
                        f'loss rate: {tree_loss_rate}, '
                        f'alpha: {tree_alpha}, '
                        f'provided_lh_fct: {tree_lh_fct} ...')
            tree = construct_tree(ls_array_names, ls_arrays, crispr_group.name,
                                  logger=logger, distance_fct=tree_distance_function,
                                  gain_rate=tree_gain_rate, loss_rate=tree_loss_rate, alpha=tree_alpha,
                                  provided_lh_fct=give_lh_fct,
                                  tree_save_path=None,
                                  tree_construction_method=tree_construction_method)
        new_tree = AdvancedTree(tree, True, model_name=crispr_group.name)

        new_dict_trees[crispr_group.name] = tree.format('newick')
        # print(ls_array_names)
        if core_genome_trees:
            # Bio.Phylo.draw_ascii(tree)
            tree, crispr_group = tree_handling(tree, crispr_group, name_inner_nodes=True, bl_eps=0)
            ls_arrays = [crispr_array.spacer_array for crispr_array in crispr_group.crispr_dict.values()]
            ls_array_names = list(crispr_group.crispr_dict.keys())
            if len(ls_arrays) < minimum_nb_of_arrays:
                logger.warning('Skipped because too few arrays after pruning.\n')
                ls_skipped_protocol.append({'name': crispr_group.name, 'repeat': crispr_group.repeat,
                                            'reason': 'Skipped because too few arrays after pruning.'})
                continue
        else:
            tree = name_inner_nodes_if_unnamed(tree)
        if len(ls_arrays) < minimum_nb_of_arrays:
            logger.warning(f'{crispr_group.name} was skipped because there is only one array.')
            ls_skipped_protocol.append({'name': crispr_group.name, 'repeat': crispr_group.repeat,
                                        'reason': 'Skipped because there is only one array.'})
            continue

        crispr_group.set_tree(tree)
        # To prevent branches from having 0 branch length
        if extend_branches:
            for node in crispr_group.tree.find_clades():
                if node.up is None:
                    node.branch_length = 0
                    continue
                if node.branch_length is None:
                    node.branch_length = 0
                node.branch_length += extend_branches
        aligned = model_tools.create_df_alignment(ls_arrays, ls_array_names)
        if core_genome_trees:
            ls_drop_cols = []
            for j, col in enumerate(aligned.columns):
                if all([c == '-' for c in aligned[col]]):
                    ls_drop_cols.append(col)
            aligned = aligned.drop(columns=ls_drop_cols)

        if provided_aligned_arrays is not None and provided_dict_duplicated_spacers is not None:
            logger.info('Importing aligned arrays...')
            resolved_df = provided_aligned_arrays[crispr_group.name]
            dict_cols_to_join = None
            mismatch_dup_stats = None
            cols = [resolved_df[col].to_list() for col in resolved_df.columns]
            uniques = [list(set(col)) for col in cols]
            raw_top_order = []
            for u in uniques:
                if len(u) > 2:
                    Exception('More than two unique values in column!')
                else:
                    raw_top_order.append(u[0])
        else:
            logger.info('Fixing mismatches (relabeling)...')
            resolved_df, raw_top_order, \
                dict_cols_to_join, mismatch_dup_stats = model_tools.get_fix_mismatches_df_alignment(aligned,
                                                                                                    )

        if sim_as_rec:
            gain_dict, loss_dict = crispr_group.gain_loss_dicts
            if hide_unobserved_spacers:
                resolved_df, raw_top_order, gain_dict, loss_dict = model_tools.remove_unobserved_spacer(resolved_df,
                                                                                                        tree,
                                                                                                        gain_dict,
                                                                                                        loss_dict,
                                                                                                        logger=logger)
            else:
                # This cleans same branch gains + losses which do not turn up in node spacers.
                for key, gains in gain_dict.items():
                    gain_dict[key] = [s for s in gains if s in raw_top_order]
                for key, losses in loss_dict.items():
                    loss_dict[key] = [s for s in losses if s in raw_top_order]
            gain_loss_dicts = (gain_dict, loss_dict)
        else:
            if hide_unobserved_spacers:
                ls_drop_cols = []
                for j, col in enumerate(resolved_df.columns):
                    if all([c == '-' for c in resolved_df[col]]):
                        ls_drop_cols.append(col)
                resolved_df = resolved_df.drop(columns=ls_drop_cols)
            gain_loss_dicts = None

        rec_m = ReconstructionTree(rec_parameter_dict, save_path, tree, True, model_name=crispr_group.name,
                                   lh_fct=give_lh_fct, sim_as_rec=sim_as_rec,
                                   sim_gain_loss_dicts=gain_loss_dicts,
                                   save_reconstructed_events=save_reconstructed_events,
                                   logger=logger,)

        tree_height = rec_m.distance_to_leafs()
        tree_length = rec_m.total_branch_length()
        min_max_tree_bl = rec_m.min_max_tree_bl()
        logger.info('Tree parameters; Distance to leaf: %s ; tree length: %s' % (tree_height,
                                                                                 tree_length))

        logger.info('Importing real data into model (also determine orders)...')
        df = rec_m.prepare_place_aligned_data_on_tree(resolved_df, raw_top_order, plot_order=False,
                                                      only_terminals=not sim_as_rec, cols_to_join=dict_cols_to_join,
                                                      dict_duplicated_spacers=provided_dict_duplicated_spacers[
                                                          crispr_group.name]
                                                      if provided_dict_duplicated_spacers is not None else None)
        avg_array_len = rec_m.avg_length_leaf_spacers()
        dict_duplicated_spacers = rec_m.get_dict_duplicated_spacers()
        # To guarantee that gain rate is very small compared to loss rate.
        rec_m.rec_model.modify_rates(rec_m.rec_model.gain_rate, rec_m.rec_model.loss_rate / tree_height,
                                     new_k=df.shape[1] * tree_height)
        # running guide reconstruction
        if not sim_as_rec:
            rec_m.ml_anc_joint(finite_bl_at_root=finite_bl_at_root)

        guide_rec_gain_dict = copy.deepcopy(rec_m.rec_gain_dict)
        nb_guide_rec_gains = sum([len(val) for val in guide_rec_gain_dict.values()])
        rec_duplications_rearrangement_placements = rec_m.duplication_rearrangement_placement(guide_rec_gain_dict)
        count_duplications_rearrangement_candidates = sum([len(value)
                                                           for value in
                                                           rec_duplications_rearrangement_placements.values()])

        guide_contra_dict = rec_m.get_rec_contradictions(find_contradictions=True)
        nb_guide_rec_contra = sum([len(val) for val in guide_contra_dict.values()])

        if plot_tree:
            dict_bg_colors = rec_m.visualize_tree(name=crispr_group.name + '_guide_rec', do_show=do_show,
                                                  path=save_path,
                                                  # sort_top_order=['9', '8', '7', '6', '5', '4', '3A', '3B', '2', '1'],
                                                  # spacer_labels_num=False,
                                                  sort_top_order='top_order',
                                                  provided_numbering=None if provided_numbering is None
                                                                             or not spacer_labels_num
                                                  else provided_numbering[crispr_group.name][0],
                                                  provided_bg_colors=None if provided_numbering is None
                                                  else provided_numbering[crispr_group.name][1],
                                                  spacer_labels_num=spacer_labels_num,
                                                  re_plot_order=False,
                                                  dpi=dpi,
                                                  figsize=figsize_rec,
                                                  )

        if not sim_as_rec:
            rec_m.ml_anc_joint_correct_contra(repeated=True, finite_bl_at_root=finite_bl_at_root)
        after_correct_contra_dict = rec_m.get_rec_contradictions(find_contradictions=True)

        count_moved_spacers = rec_m.count_moved_spacers(guide_rec_gain_dict)

        # Categorize duplications and rearrangements.
        rec_duplications, rec_rearrangements, rec_reacquisitions, \
            rec_default_duplications, rec_other_dup_events = rec_m.distinction_dup_rearrangement()
        count_duplications, count_rearrangements, \
            count_reacquisitions, count_ind_dups_default, \
            count_other_dup_events, \
            count_indep_acquisitions_not_dup_candidates = rec_m.get_rec_duplication_type_counts()

        logger.info(f'Length of topological order: {len(rec_m.top_order)}')
        nb_gained_spacers = len(rec_m.top_order)

        nb_rec_insertions = sum([len(val) for val in rec_m.rec_gain_dict.values()])
        nb_rec_deletions = sum([len(val) for val in rec_m.rec_loss_dict.values()])
        nb_unique_spacers = len(rec_m.top_order) - count_duplications_rearrangement_candidates
        dict_matrices = rec_m.generate_possible_losses_sets()
        dict_max_length_losses, ls_all_max_lengths, ls_all_max_lengths_norm = rec_m.get_max_length_loss_sets()
        ls_rel_loss_pos, ls_losses, \
            ls_nb_ex_spacers, avg_nb_ex_spacers, \
            (fes_m_presence, les_m_presence, l_s_m_presence, f_s_m_presence), \
            (global_fes_m_presence, global_les_m_presence, global_l_s_presence) = rec_m.get_relative_loss_stats()

        if give_lh_fct is not None:
            lh_0, lh_1, opt_result_0, opt_result_1 = rec_m.compute_lh_ratio_for_given_lh_fct()
        else:
            lh_0, lh_1, opt_result_0, opt_result_1 = rec_m.compute_lh_ratio(alpha_bias_correction=alpha_bias_correction,
                                                                            rho_bias_correction=rho_bias_correction, )

        lh_ratio = lh_0 / lh_1 if lh_1 != 0 else np.nan
        ln_lh_0, ln_lh_1 = -opt_result_0.fun, -opt_result_1.fun
        ln_lh_ratio = 2 * (ln_lh_1 - ln_lh_0)
        test_result, quantile = stats.test_significance_ratio_chi2(ln_lh_ratio, significance_level)
        preferred_model = 'BDM' if test_result else 'IDM'
        if nb_rec_deletions == 0:
            logger.info('No deletions were reconstructed, therefore the test results are not meaningful '
                        '(especially estimated parameters).')
        logger.info(f'ln(lh_idm): {ln_lh_0} ; ln(lh_bdm): {ln_lh_1} ; '
                    f'test statistic (2*(ln(lh_bdm) - ln(lh_idm)): {ln_lh_ratio} ; '
                    f'test significant?: {test_result} ; '
                    f'chi2 quantile: {quantile} ; '
                    f'of chosen significance level: {significance_level}'
                    )
        logger.info(f'Test decides for: {preferred_model} ; '
                    f'estimated ML IDM parameter: {opt_result_0.x} ; '
                    f'estimated ML BDM parameters [deletion_rate, alpha]: {opt_result_1.x}')

        # alternative estimators with filtered joined losses
        if not alternative_parameter_estimation or len(rec_m.rec_gain_dict[rec_m.root.name]) == 0:
            ex_y_s_res_lr_0, ex_y_s_res_lr_1, ex_y_s_res_alpha_1 = np.nan, np.nan, np.nan
            con_o_s_res_lr_0, con_o_s_res_lr_1, con_o_s_res_alpha_1 = np.nan, np.nan, np.nan
            if len(rec_m.rec_gain_dict[rec_m.root.name]) == 0:
                logger.warning('The root has NO insertions! Therefore alternative parameter estimations based on root '
                               'spacers were set to np.nan.')
        else:
            if give_lh_fct is not None:
                _, _, ex_y_s_res_0, ex_y_s_res_1 = rec_m.compute_lh_ratio_for_given_lh_fct(
                    filter_by='exclude_young_spacers')
                _, _, con_o_s_res_0, con_o_s_res_1 = rec_m.compute_lh_ratio_for_given_lh_fct(
                    filter_by='contains_old_spacers')
            else:
                _, _, ex_y_s_res_0, ex_y_s_res_1 = rec_m.compute_lh_ratio(filter_by='exclude_young_spacers',
                                                                          alpha_bias_correction=alpha_bias_correction,
                                                                          rho_bias_correction=rho_bias_correction)
                _, _, con_o_s_res_0, con_o_s_res_1 = rec_m.compute_lh_ratio(filter_by='contains_old_spacers',
                                                                            alpha_bias_correction=alpha_bias_correction,
                                                                            rho_bias_correction=rho_bias_correction)
            ex_y_s_res_lr_0, ex_y_s_res_lr_1 = ex_y_s_res_0.x, ex_y_s_res_1.x[0]
            ex_y_s_res_alpha_1 = ex_y_s_res_1.x[1]
            con_o_s_res_lr_0, con_o_s_res_lr_1 = con_o_s_res_0.x, con_o_s_res_1.x[0]
            con_o_s_res_alpha_1 = con_o_s_res_1.x[1]

        if plot_tree:
            dict_bg_colors = rec_m.visualize_tree(name=crispr_group.name + '_rec', do_show=False,
                                                  path=save_path, sort_top_order='nb_leafs', indicate_joint_del=True,
                                                  provided_numbering=None if provided_numbering is None
                                                                             or not spacer_labels_num
                                                  else provided_numbering[crispr_group.name][0],
                                                  provided_bg_colors=None if provided_numbering is None
                                                  else provided_numbering[crispr_group.name][1],
                                                  spacer_labels_num=spacer_labels_num,
                                                  re_plot_order=False,
                                                  dpi=dpi,
                                                  figsize=figsize_rec,
                                                  )
        else:
            dict_bg_colors = None
        if plot_order:
            rec_m.visualize_order(save_path, crispr_group.name,
                                  save_path_dot=os.path.join(crispr_group.name + '_order.dot'),
                                  color_dict=dict_bg_colors,
                                  spacer_or_number='number',
                                  provided_numbering=None if provided_numbering is None
                                  else provided_numbering[crispr_group.name][0]
                                  )
        dict_array_lengths = rec_m.get_leaf_array_lengths()

        run_time = run_timer.time_from_last_checkpoint()

        dict_protocol = {'name': crispr_group.name,
                         'deletion_rate_bdm': opt_result_1.x[0], 'alpha_bdm': opt_result_1.x[1],
                         'per_spacer_deletion_rate_bdm': opt_result_1.x[0] * opt_result_1.x[1],
                         'insertion_rate_bdm': opt_result_1.x[0] * opt_result_1.x[1] * avg_array_len,
                         'deletion_rate_idm': opt_result_0.x,
                         'insertion_rate_idm': opt_result_0.x * avg_array_len,
                         'ln_lh_idm': -opt_result_0.fun, 'ln_lh_bdm': -opt_result_1.fun,
                         'test_statistic (-2*ln_lh_ratio)': ln_lh_ratio, 'test result': test_result,
                         'Deletion model preferred by LRT': preferred_model,
                         'chi2_quantile': quantile,
                         'significance level': significance_level,
                         'avg array length': avg_array_len, 'nb of spacers in alignment': nb_gained_spacers,
                         'nb of unique spacers': nb_unique_spacers,
                         'nb of reconstructed insertions': nb_rec_insertions,
                         'nb of reconstructed deletions': nb_rec_deletions,
                         'nb of reconstructed dup/rearr candidates': count_duplications_rearrangement_candidates,
                         'nb of reconstructed rearrangements': count_rearrangements,
                         'nb of reconstructed duplications': count_duplications,
                         'nb of reconstructed reacquisitions': count_reacquisitions,
                         'nb of reconstructed independent gains': count_ind_dups_default,
                         'nb of reconstructed other dup. events': count_other_dup_events,
                         'nb of reconstructed ind. acquisitions, not duplicates': count_indep_acquisitions_not_dup_candidates,
                         'tree height': tree_height,
                         'tree length': tree_length,
                         'min/max tree branch lengths': min_max_tree_bl,
                         'nb of leafs (before combining non-uniques)': prev_len,
                         'nb of leafs (after combining non-uniques)': len(ls_arrays),
                         'nb of unique spacer arrays': nb_unique_spacer_arrays,
                         'nb of guide reconstruction insertions': nb_guide_rec_gains,
                         'nb of guide reconstruction contradictions': nb_guide_rec_contra,
                         'nb of moved spacers by refinement': count_moved_spacers,
                         'run_time': run_time,
                         'array names': ls_array_names, 'combined array names': dict_combined_arrays,
                         'leaf array lengths': dict_array_lengths,
                         'lh_idm': lh_0, 'lh_bdm': lh_1,
                         'excluding young spacers deletion_rate_idm': ex_y_s_res_lr_0,
                         'excluding young spacers deletion_rate_bdm': ex_y_s_res_lr_1,
                         'excluding young spacers alpha_bdm': ex_y_s_res_alpha_1,
                         'containing old spacers deletion_rate_idm': con_o_s_res_lr_0,
                         'containing old spacers deletion_rate_bdm': con_o_s_res_lr_1,
                         'containing old spacers alpha_bdm': con_o_s_res_alpha_1,
                         'relative deletion positions': ls_rel_loss_pos,
                         'nb of existent spacers': ls_nb_ex_spacers,
                         'all max block deletion lengths': ls_all_max_lengths,
                         'all max block deletion lengths (normalized)': ls_all_max_lengths_norm,
                         'fes +m presence': fes_m_presence, 'les -m presence': les_m_presence,
                         'ls -m presence': l_s_m_presence,
                         'fs +m presence': f_s_m_presence,
                         'global fes +m presence': global_fes_m_presence,
                         'global les -m presence': global_les_m_presence,
                         'global ls -m presence': global_l_s_presence,
                         }
        if metadata:
            dict_crispr_group_metadata = {'repeat': crispr_group.repeat,
                                          'crispr type': crispr_group.cas_type,
                                          'chromosome_plasmid': crispr_group.chromosome_plasmid,
                                          'group species': crispr_group.group_species,
                                          'species by array': [a.species.species
                                                               for a in crispr_group.crispr_dict.values()
                                                               if a.species is not None],
                                          'genus by array': [a.species.genus
                                                             for a in crispr_group.crispr_dict.values()
                                                             if a.species is not None],
                                          'kingdom': crispr_group.kingdom,
                                          'array orientation': [a.orientation
                                                                for a in crispr_group.crispr_dict.values()],
                                          }
            dict_protocol.update(dict_crispr_group_metadata)
        if filter_fct(dict_protocol):
            ls_dict_protocol.append(dict_protocol)

            dict_data_for_lh_ratio[crispr_group.name] = rec_m.get_data_for_lh_ratio(filter_by=False)
            if i % 10 == 1:
                df_protocol = pd.DataFrame(ls_dict_protocol)
                df_protocol.to_csv(os.path.join(save_path, '_'.join(['0_protocol']) + '.csv'))
                df_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_protocol']) + '.pkl'))
        else:
            # for bugfixing
            # _ = rec_m.visualize_tree(name=crispr_group.name + '_boring', do_show=False,
            #                          path=save_path, sort_top_order='nb_leafs', indicate_joint_del=True,
            #                          provided_numbering=None)
            logger.info(f'{crispr_group.name} seems to be trivial data.')
            ls_boring_protocol.append(dict_protocol)

        new_dict_provided_numberings[crispr_group.name] = (rec_m.spacer_names_to_numbers, dict_bg_colors)
        new_dict_provided_aligned_arrays[crispr_group.name] = df
        new_dict_provided_duplicated_spacers[crispr_group.name] = dict_duplicated_spacers
        logger.info(f'Run time of group "{crispr_group.name}": {run_time} s\n')
    if tree_save_path is not None:
        with open(os.path.join(tree_save_path, 'dict_nwk_trees.pkl'), 'wb') as f:
            pickle.dump(new_dict_trees, f)

    df_protocol = pd.DataFrame(ls_dict_protocol)
    df_protocol.to_csv(os.path.join(save_path, '_'.join(['0_protocol']) + '.csv'))
    df_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_protocol']) + '.pkl'))
    protocol_file_name = os.path.join(save_path, '_'.join(['0_protocol']))

    df_boring_protocol = pd.DataFrame(ls_boring_protocol)
    df_boring_protocol.to_csv(os.path.join(save_path, '_'.join(['0_protocol_wo_trivial']) + '.csv'))
    df_boring_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_protocol_wo_trivial']) + '.pkl'))

    df_skipped_protocol = pd.DataFrame(ls_skipped_protocol)
    df_skipped_protocol.to_csv(os.path.join(save_path, '_'.join(['0_protocol_skipped']) + '.csv'))
    df_skipped_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_protocol_skipped']) + '.pkl'))

    # with open(os.path.join(save_path, '0_rec_m_saved.pkl'), 'rb') as f:
    #     ls_rec_m = import_data.unpickle(f)
    # ls_rec_m.sort(key=lambda x: x[0])
    # ls_rec_m = [y[1] for y in ls_rec_m]
    if group_by is not None:
        logger.info(f'Grouping trees for parameter estimation by: {group_by}')
        if isinstance(group_by, int):
            val_group_by = []
            for i in range(df_protocol.shape[0] // group_by):
                val_group_by += [i] * group_by
            val_group_by += [max(val_group_by)] * (df_protocol.shape[0] % group_by)
        elif isinstance(group_by, pd.Series):
            val_group_by = [group_by[name] for name in df_protocol['name']]
        else:
            val_group_by = df_protocol[group_by]
        ls_data_for_lh_ratio = [dict_data_for_lh_ratio[name] for name in df_protocol['name']]
        group_names, groups_lh_0, groups_lh_1, \
            groups_result_0, groups_result_1 = pooling_for_parameter_estimation(ls_data_for_lh_ratio,
                                                                                group_by=val_group_by,
                                                                                give_lh_fct=give_lh_fct)

        ls_loss_rates_0 = [g.x for g in groups_result_0]
        ls_loss_rates_1 = [g.x[0] for g in groups_result_1]
        ls_alphas_1 = [g.x[1] for g in groups_result_1]
        ls_ln_lh_ratio = [2 * (g_0.fun - g_1.fun) for g_0, g_1 in zip(groups_result_0, groups_result_1)]
        ls_test_results = []
        ls_confidences = []
        ls_sig_values = []
        for llh in ls_ln_lh_ratio:
            test_result, confidence = stats.test_significance_ratio_chi2(llh, significance_level)
            ls_test_results.append(test_result)
            ls_confidences.append(confidence)
            ls_sig_values.append(significance_level)

        df_group_protocol = pd.DataFrame(list(zip(group_names, groups_lh_0, groups_lh_1,
                                                  ls_loss_rates_0, ls_loss_rates_1,
                                                  ls_alphas_1, ls_ln_lh_ratio, ls_test_results,
                                                  ls_confidences, ls_sig_values,
                                                  )),
                                         columns=['group_by', 'lh_0', 'lh_1', 'loss_rate_0', 'loss_rate_1',
                                                  'alpha_1', '-2*ln_lh_ratio', 'test result', 'test_confidence',
                                                  'significance value'])
        df_group_protocol.to_csv(os.path.join(save_path, '_'.join(['0_group_by_protocol']) + '.csv'))
        df_group_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_group_by_protocol']) + '.pkl'))

    # Should be included as group by case.
    if len(dict_crispr_groups) > 1:
        logger.info('Computing likelihood ratio and parameter estimates over all groups...')
        lh_0, lh_1, opt_result_0, opt_result_1 = model_tools.compute_lh_ratio_of_multiple_trees(
            list(dict_data_for_lh_ratio.values()),
            method='Nelder-Mead',
            filter_by=False,
            lambdifyed_lh_fct=give_lh_fct)

        ln_lh_0, ln_lh_1 = -opt_result_0.fun, -opt_result_1.fun
        ln_lh_ratio = 2 * (ln_lh_1 - ln_lh_0)
        test_result, quantile = stats.test_significance_ratio_chi2(ln_lh_ratio, significance_level)
        preferred_model = 'BDM' if test_result else 'IDM'
        logger.info(f'LRT and parameter estimates overall groups: \n'
                    f'ln(lh_idm): {ln_lh_0} ; ln(lh_bdm): {ln_lh_1} ; '
                    f'test statistic (2*(ln(lh_bdm) - ln(lh_idm)): {ln_lh_ratio} ; '
                    f'test significant?: {test_result} ; '
                    f'chi2 quantile: {quantile} ; '
                    f'of chosen significance level: {significance_level}')
        logger.info(f'Test decides for: {preferred_model} ; '
                    f'estimated ML IDM parameter: {opt_result_0.x} ; '
                    f'estimated ML BDM parameters [deletion_rate, alpha]: {opt_result_1.x}')

    logger.info(f'Runtime of all experiments: {run_timer.time_from_start():.2f} s ; detailed results saved to '
                f'{protocol_file_name}')

    if dict_trees and use_provided_trees:
        new_dict_trees = dict_trees

    return df_protocol, df_boring_protocol, (dict_data_for_lh_ratio, give_lh_fct), new_dict_provided_numberings, \
        new_dict_trees, new_dict_provided_aligned_arrays, new_dict_provided_duplicated_spacers