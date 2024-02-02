import logging
import os
import pandas as pd
import copy
import pickle
import numpy as np
import random
import json

from model import mafft_integration, distance_tree_constructor, model_tools
from model.helpers import misc, import_data, stats
from model.orientation import orientation_tools
from model.bdm_likelihood_computation import symbolic_lh_computation
from model.experiments_tools import expand_sim_parameters, multiple_selection_fcts, construct_tree, filter_fct, \
    load_align_omer_data, load_arrays_from_pkl, align_crispr_groups_for_sim_as_rec, pooling_for_parameter_estimation, \
    load_align_single_fasta, expand_sim_parameters_based_on_import, tree_handling
from model.simulation_tree import SimulationTree
from model.reconstruction_tree import ReconstructionTree
from model.data_classes.advanced_tree import AdvancedTree
from model.summary import write_summary, compose_summary_dict

WORK_PATH = os.path.join('data', 'simulation_alignment')
TREE_GENERATION_PARAMETERS = {
    'simple_insertion_rate': 6522.48,
    'simple_deletion_rate': 137.20,
    'simple_alpha': 2.7278,
    'precise_insertion_rate': 0.6108,
    'precise_deletion_rate': 0.1830,
    'precise_alpha': 3.3377,
}
LH_FCT_PATH = os.path.join('data', '0_lh', 'death_lh_up_to_100_lambdified.pkl')


def run_multiple_groups(ls_data_path, save_path, rec_parameter_dict, lh_fct=None, logger=None, plot_tree=True,
                        do_show=False, combine_non_unique_arrays=False, determine_orientation=True,
                        orientation_decision_boundary=10,
                        tree_path=None, plot_order=True, significance_level=0.05, extend_branches=False,
                        tree_lh_fct=None, tree_insertion_rate=None, tree_deletion_rate=None, tree_alpha=None,
                        alpha_bias_correction=True, rho_bias_correction=True,
                        ):
    dict_crispr_groups = {}
    run_timer = misc.RunTimer()
    save_data_folder = os.path.join(save_path, 'additional_data')
    logger.info(f'Loading and aligning the data...')
    for i, data_path in enumerate(ls_data_path):
        group_name = os.path.split(data_path)[-1].split('.')[0]
        ind_save_path = os.path.join(save_data_folder, 'work_folder', group_name)
        if not os.path.exists(ind_save_path):
            os.makedirs(ind_save_path)
        logger.info(f'{i + 1} / {len(ls_data_path)} aligning group {group_name}... \n')
        crispr_group = load_align_single_fasta(data_path, ind_save_path,
                                               group_name=group_name,
                                               logger=logger,
                                               mafft_options=None,
                                               save_path=None)
        dict_crispr_groups[group_name] = crispr_group
    # logger.debug(f'Loading and aligning time: {run_timer.time_from_last_checkpoint()}')
    print_save_path = os.path.join(save_data_folder, 'work_folder')
    logger.debug(f'Saving CRISPR groups to: {print_save_path}')
    if not os.path.exists(save_data_folder):
        os.makedirs(save_data_folder)
    pickle.dump(dict_crispr_groups, open(os.path.join(save_data_folder, 'dict_crispr_groups.pkl'), 'wb'))

    if lh_fct.lower() in ['ode_based', 'precise']:
        logger.info(f'Loading lh function from {LH_FCT_PATH}')
        # Implement loading a shortened likelihood function depending on something.
        lh_fct = symbolic_lh_computation.load_pickled_lambdified_lh_fct(LH_FCT_PATH)
        logger.info(f'Loading time: {run_timer.time_from_last_checkpoint()}')

    if tree_lh_fct in ['ode_based', 'precise']:
        tree_lh_fct = symbolic_lh_computation.load_pickled_lambdified_lh_fct(LH_FCT_PATH) if lh_fct not in ['ode_based',
                                                                                                            'precise'] \
            else lh_fct
    dict_crispr_groups_for_reverse = copy.deepcopy(dict_crispr_groups) if determine_orientation else {}

    logger.info(f'Beginning reconstruction of (forward oriented) group(s)...\n')
    df_rec_protocol, \
        df_rec_protocol_trivial, _, \
        dict_provided_numbering, dict_trees_forward, \
        dict_provided_aligned_arrays, \
        dict_provided_duplicated_spacers = run_reconstruction(rec_parameter_dict,
                                                              dict_crispr_groups,
                                                              save_path=os.path.join(
                                                                  save_path,
                                                                  '0_forward'),
                                                              plot_tree=plot_tree,
                                                              lh_fct=lh_fct,
                                                              logfile_path=None,
                                                              logger=logger, do_show=do_show,
                                                              combine_non_unique_arrays=combine_non_unique_arrays,
                                                              tree_path=tree_path,
                                                              use_provided_trees=True if tree_path is not None
                                                              else False,
                                                              tree_save_path=None,
                                                              hide_unobserved_spacers=False,
                                                              plot_order=plot_order,
                                                              finite_bl_at_root=True,
                                                              group_by=None,
                                                              significance_level=significance_level,
                                                              extend_branches=extend_branches,
                                                              tree_distance_function=tree_lh_fct,
                                                              tree_gain_rate=tree_insertion_rate,
                                                              tree_loss_rate=tree_deletion_rate,
                                                              tree_alpha=tree_alpha,
                                                              alpha_bias_correction=alpha_bias_correction,
                                                              rho_bias_correction=rho_bias_correction,
                                                              core_genome_trees=False,
                                                              metadata=False,
                                                              )
    dict_trees = dict_trees_forward
    df_rec_protocol = df_rec_protocol.set_index('name')
    if not df_rec_protocol_trivial.empty:
        df_rec_protocol_trivial = df_rec_protocol_trivial.set_index('name')

    for gn, provided_numbering in dict_provided_numbering.items():
        dict_spacer_names_to_numbers = provided_numbering[0]
        p = os.path.join(save_data_folder, 'spacer_fasta', gn + '_spacer_name_to_sp_number.json')
        if not os.path.exists(os.path.dirname(p)):
            os.makedirs(os.path.dirname(p))
        json.dump(dict_spacer_names_to_numbers, open(p,
                                                     'w'))
    if determine_orientation:
        logger.info(f'Beginning reconstruction of reversed group(s)...\n')
        if not os.path.exists(os.path.join(save_path, '0_reversed')):
            os.makedirs(os.path.join(save_path, '0_reversed'))

        for gn, crispr_group in dict_crispr_groups_for_reverse.items():
            for n, ca in crispr_group.crispr_dict.items():
                ca.spacer_array = ca.spacer_array[::-1]
        dict_provided_aligned_arrays = {k: v[v.columns[::-1]] for k, v in dict_provided_aligned_arrays.items()}

        df_rec_protocol_reversed, \
            df_rec_protocol_reversed_trivial, \
            _, _, dict_trees_reversed, \
            _, _ = run_reconstruction(rec_parameter_dict,
                                      dict_crispr_groups_for_reverse,
                                      save_path=os.path.join(
                                          save_path,
                                          '0_reversed'),
                                      plot_tree=plot_tree,
                                      lh_fct=lh_fct,
                                      logfile_path=None,
                                      logger=logger, do_show=do_show,
                                      combine_non_unique_arrays=combine_non_unique_arrays,
                                      tree_path=tree_path,
                                      use_provided_trees=True if tree_path is not None
                                      else False,
                                      tree_save_path=None,
                                      hide_unobserved_spacers=False,
                                      plot_order=plot_order,
                                      finite_bl_at_root=True,
                                      group_by=None,
                                      significance_level=significance_level,
                                      extend_branches=extend_branches,
                                      tree_distance_function=tree_lh_fct,
                                      tree_gain_rate=tree_insertion_rate,
                                      tree_loss_rate=tree_deletion_rate,
                                      tree_alpha=tree_alpha,
                                      alpha_bias_correction=alpha_bias_correction,
                                      rho_bias_correction=rho_bias_correction,
                                      core_genome_trees=False,
                                      provided_numbering=dict_provided_numbering,
                                      provided_aligned_arrays=dict_provided_aligned_arrays,
                                      provided_dict_duplicated_spacers=dict_provided_duplicated_spacers,
                                      metadata=False)
        dict_trees = {'forward': dict_trees_forward,
                      'reversed': dict_trees_reversed}

        df_rec_protocol_reversed = df_rec_protocol_reversed.set_index('name')
        if not df_rec_protocol_reversed_trivial.empty:
            df_rec_protocol_reversed_trivial = df_rec_protocol_reversed_trivial.set_index('name')

        df_rec_protocol, df_oriented_rec_protocol, dict_trees = orientation_tools.compare_likelihoods_for_orientation(
            df_rec_protocol,
            df_rec_protocol_reversed,
            df_rec_protocol_trivial,
            df_rec_protocol_reversed_trivial,
            decision_boundary=orientation_decision_boundary, dict_trees=dict_trees)

        detailed_protocols_save_path = os.path.join(save_path, 'detailed_results')
        if not os.path.exists(detailed_protocols_save_path):
            os.makedirs(detailed_protocols_save_path)
        df_rec_protocol.to_csv(os.path.join(detailed_protocols_save_path,
                                            '_'.join(['0_detailed_results_w_orientation']) + '.csv'))
        df_rec_protocol.to_pickle(os.path.join(detailed_protocols_save_path,
                                               '_'.join(['0_detailed_results_w_orientation']) + '.pkl'))
        df_oriented_rec_protocol.to_csv(os.path.join(detailed_protocols_save_path,
                                                     '0_detailed_results_wo_trivial_groups.csv'))
        df_oriented_rec_protocol.to_pickle(os.path.join(detailed_protocols_save_path,
                                                        '0_detailed_results_wo_trivial_groups.pkl'))
    df_results_wo_details = df_rec_protocol.drop(columns=['relative deletion positions', 'nb of existent spacers',
                                                          'all max block deletion lengths',
                                                          'all max block deletion lengths (normalized)',
                                                          'reversed_lh_idm', 'reversed_lh_bdm',
                                                          'lh_idm / reversed_lh_idm',
                                                          'lh_bdm / reversed_lh_bdm',
                                                          'lh_idm', 'lh_bdm',
                                                          'excluding young spacers deletion_rate_idm',
                                                          'excluding young spacers deletion_rate_bdm',
                                                          'excluding young spacers alpha_bdm',
                                                          'containing old spacers deletion_rate_idm',
                                                          'containing old spacers deletion_rate_bdm',
                                                          'containing old spacers alpha_bdm', 'fes +m presence',
                                                          'les -m presence',
                                                          'ls -m presence', 'fs +m presence', 'global fes +m presence',
                                                          'global les -m presence', 'global ls -m presence',
                                                          'combined array names',
                                                          ],
                                                 errors='ignore')

    df_results_wo_details.to_csv(os.path.join(save_path, '_'.join(['0_results']) + '.csv'))

    # write_summary()
    # summary = {**summary_forward, **summary_reversed}

    tree_save_path = os.path.join(save_path, 'additional_data')
    if tree_save_path is not None:
        with open(os.path.join(tree_save_path, 'dict_nwk_trees.json'), 'w') as f:
            json.dump(dict_trees, f)

    return df_results_wo_details


def run_pickled_data(rec_parameter_dict, data_path, save_path=None, plot_tree=True,
                     logfile_path=None,
                     do_show=False,
                     combine_non_unique_arrays=True,
                     tree_path=None,
                     tree_save_path=None,
                     hide_unobserved_spacers=False,
                     selection_fct=None,
                     plot_order=True,
                     lh_fct=None,
                     significance_level=0.05,
                     finite_bl_at_root=True,
                     group_by=None,
                     alignment_options_dict=None,
                     extend_branches=False,
                     tree_distance_function='breakpoint',
                     tree_gain_rate=None,
                     tree_loss_rate=None,
                     tree_alpha=None,
                     determine_orientation=False,
                     orient_boundary=0,
                     alpha_bias_correction=False,
                     rho_bias_correction=False,
                     core_genome_trees=False,
                     ):
    if not os.path.exists(os.path.split(logfile_path)[0]):
        os.makedirs(os.path.split(logfile_path)[0])

    if isinstance(logfile_path, str):
        logger = misc.create_logger('pickled data based reconstruction', logging.INFO, outfile=logfile_path)
    else:
        logger = logfile_path

    if data_path.split('.')[-1] in ['pkl', 'pickle']:
        logger.info(f'Loading data from pickle {data_path}')
        dict_crispr_groups = load_arrays_from_pkl(data_path)
        logger.info(f'Number of groups before reconstruction: {len(dict_crispr_groups)}')
        nb_arrays = sum([len(g.crispr_dict) for g in dict_crispr_groups.values()])
        logger.info(f'Number of arrays before reconstruction: {nb_arrays}')
    else:
        logger.info(f'Loading and aligning from files {data_path}')
        dict_crispr_groups = load_align_omer_data(data_path, mafft_options=alignment_options_dict['mafft_options'],
                                                  exclude_files=alignment_options_dict['exclude_files'],
                                                  save_path=alignment_options_dict['save_path'])

    dict_crispr_groups_for_reverse = copy.deepcopy(dict_crispr_groups) if determine_orientation else {}
    (df_rec_protocol, df_rec_protocol_boring, (dict_data_for_group_by, lh_fct), dict_provided_numbering,
     dict_trees_forward,
     dict_provided_aligned_arrays,
     dict_provided_duplicated_spacers) = run_reconstruction(rec_parameter_dict,
                                                            dict_crispr_groups,
                                                            save_path=os.path.join(
                                                                save_path,
                                                                '0_forward'),
                                                            plot_tree=plot_tree,
                                                            lh_fct=lh_fct,
                                                            logfile_path=logfile_path,
                                                            logger=logger, do_show=do_show,
                                                            combine_non_unique_arrays=combine_non_unique_arrays,
                                                            tree_path=tree_path,
                                                            use_provided_trees=True if tree_path is not None else False,
                                                            tree_save_path=os.path.join(
                                                                tree_save_path,
                                                                '0_forward'),
                                                            hide_unobserved_spacers=hide_unobserved_spacers,
                                                            selection_fct=selection_fct,
                                                            plot_order=plot_order,
                                                            finite_bl_at_root=finite_bl_at_root,
                                                            group_by=group_by,
                                                            significance_level=significance_level,
                                                            extend_branches=extend_branches,
                                                            tree_distance_function=tree_distance_function,
                                                            tree_gain_rate=tree_gain_rate,
                                                            tree_loss_rate=tree_loss_rate,
                                                            tree_alpha=tree_alpha,
                                                            alpha_bias_correction=alpha_bias_correction,
                                                            rho_bias_correction=rho_bias_correction,
                                                            core_genome_trees=core_genome_trees, )
    dict_trees = dict_trees_forward
    df_rec_protocol = df_rec_protocol.set_index('name')
    df_rec_protocol_boring = df_rec_protocol_boring.set_index('name')
    if determine_orientation:
        if not os.path.exists(os.path.join(save_path, '0_reversed')):
            os.makedirs(os.path.join(save_path, '0_reversed'))

        for gn, crispr_group in dict_crispr_groups_for_reverse.items():
            for n, ca in crispr_group.crispr_dict.items():
                ca.spacer_array = ca.spacer_array[::-1]
        dict_provided_aligned_arrays = {k: v[v.columns[::-1]] for k, v in dict_provided_aligned_arrays.items()}

        (df_rec_protocol_reversed,
         df_rec_protocol_reversed_boring,
         (dict_data_for_group_by_reversed, lh_fct), _, dict_trees_reversed,
         _, _) = run_reconstruction(rec_parameter_dict,
                                    dict_crispr_groups_for_reverse,
                                    save_path=os.path.join(
                                        save_path,
                                        '0_reversed'),
                                    plot_tree=plot_tree,
                                    lh_fct=lh_fct,
                                    logfile_path=logfile_path,
                                    logger=logger,
                                    do_show=do_show,
                                    combine_non_unique_arrays=combine_non_unique_arrays,
                                    tree_path=tree_path,
                                    use_provided_trees=True
                                    if tree_path is not None else False,
                                    tree_save_path=os.path.join(
                                        tree_save_path,
                                        '0_reversed'),
                                    hide_unobserved_spacers=hide_unobserved_spacers,
                                    selection_fct=selection_fct,
                                    plot_order=plot_order,
                                    finite_bl_at_root=finite_bl_at_root,
                                    group_by=group_by,
                                    significance_level=significance_level,
                                    extend_branches=extend_branches,
                                    tree_distance_function=tree_distance_function,
                                    tree_gain_rate=tree_gain_rate,
                                    tree_loss_rate=tree_loss_rate,
                                    tree_alpha=tree_alpha,
                                    provided_numbering=dict_provided_numbering,
                                    alpha_bias_correction=alpha_bias_correction,
                                    rho_bias_correction=rho_bias_correction,
                                    core_genome_trees=core_genome_trees,
                                    provided_aligned_arrays=dict_provided_aligned_arrays,
                                    provided_dict_duplicated_spacers=dict_provided_duplicated_spacers, )
        dict_trees = {'forward': dict_trees_forward,
                      'reversed': dict_trees_reversed}
        df_rec_protocol_reversed = df_rec_protocol_reversed.set_index('name')
        df_rec_protocol_reversed_boring = df_rec_protocol_reversed_boring.set_index('name')

        df_rec_protocol, df_oriented_rec_protocol, dict_trees = orientation_tools.compare_likelihoods_for_orientation(
            df_rec_protocol,
            df_rec_protocol_reversed,
            df_rec_protocol_boring,
            df_rec_protocol_reversed_boring,
            decision_boundary=orient_boundary, dict_trees=dict_trees)

        if group_by is not None:
            logger.info(f'Grouping trees for parameter estimation by: {group_by}')
            if isinstance(group_by, int):
                val_group_by = []
                for i in range(df_oriented_rec_protocol.shape[0] // group_by):
                    val_group_by += [i] * group_by
                val_group_by += [max(val_group_by)] * (df_oriented_rec_protocol.shape[0] % group_by)
            elif isinstance(group_by, pd.Series):
                # solve this issue here pooling and group_by do not agree in dimensions!!!!
                val_group_by = group_by
            else:
                val_group_by = df_oriented_rec_protocol[group_by]
            ls_data_for_lh_ratio = []
            for idx in df_oriented_rec_protocol.index:
                # if idx in list(df_rec_protocol_boring['name']) \
                #         and idx in list(df_rec_protocol_reversed_boring['name']):
                #     continue
                ls_data_for_lh_ratio.append(dict_data_for_group_by_reversed[idx] \
                                                if df_oriented_rec_protocol.loc[
                                                       idx, 'our predicted orientation'] == 'Reverse' \
                                                else dict_data_for_group_by[idx])
            group_names, groups_lh_0, groups_lh_1, \
                groups_result_0, groups_result_1 = pooling_for_parameter_estimation(ls_data_for_lh_ratio,
                                                                                    group_by=val_group_by,
                                                                                    give_lh_fct=lh_fct)

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
                                             columns=['group_by', 'lh_idm', 'lh_bdm', 'deletion_rate_idm',
                                                      'deletion_rate_bdm',
                                                      'alpha_bdm', 'test statistic (-2*ln_lh_ratio)',
                                                      'test result', 'chi2_quantile',
                                                      'significance level'])
            df_group_protocol.to_csv(os.path.join(save_path, '0_final_group_by_protocol.csv'))

        df_rec_protocol.to_csv(os.path.join(save_path, '_'.join(['0_protocol_w_orientation']) + '.csv'))
        df_rec_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_protocol_w_orientation']) + '.pkl'))
        df_oriented_rec_protocol.to_csv(os.path.join(save_path, '0_final_oriented_protocol.csv'))
        df_oriented_rec_protocol.to_pickle(os.path.join(save_path, '0_final_oriented_protocol.pkl'))
    if tree_save_path is not None:
        with open(os.path.join(tree_save_path, 'final_dict_nwk_trees.json'), 'w') as f:
            json.dump(dict_trees, f)
    return


def run_simulation_and_reconstruction(sim_parameter_dict, rec_parameter_dict, save_path=None,
                                      rec_save_path=None, plot_tree=True,
                                      sim_save_path=None, sim_plot_tree=True,
                                      sim_logfile_path=None,
                                      logfile_path=None,
                                      do_show=False,
                                      sim_as_rec=False,
                                      hide_unobserved_spacers=False,
                                      selection_fct=None,
                                      plot_order=True,
                                      significance_level=0.05,
                                      alignment_options=None,
                                      load_sim_from_pkl=False,
                                      finite_bl_at_root=True,
                                      sim_simple_alignment=False,
                                      group_by=None,
                                      extend_branches=False,
                                      use_sim_tree=True,
                                      tree_distance_function='breakpoint',
                                      tree_gain_rate=None,
                                      tree_loss_rate=None,
                                      tree_alpha=None,
                                      tree_save_path=None,
                                      determine_orientation=False,
                                      randomize_orientation=False,
                                      orient_decision_boundary=0,
                                      load_sim_param_from_pkl=None,
                                      alpha_bias_correction=False,
                                      rho_bias_correction=False,
                                      ):
    """
    This function runs a simulation and a reconstruction and compares the results.
    :param rho_bias_correction:
    :param alpha_bias_correction:
    :param sim_parameter_dict:
    :param rec_parameter_dict:
    :param save_path:
    :param rec_save_path:
    :param plot_tree:
    :param sim_save_path:
    :param sim_plot_tree:
    :param sim_logfile_path:
    :param logfile_path:
    :param do_show:
    :param sim_as_rec:
    :param hide_unobserved_spacers:
    :param selection_fct:
    :param plot_order:
    :param significance_level:
    :param alignment_options:
    :param load_sim_from_pkl:
    :param finite_bl_at_root:
    :param sim_simple_alignment:
    :param group_by:
    :param extend_branches:
    :param use_sim_tree:
    :param tree_distance_function:
    :param tree_gain_rate:
    :param tree_loss_rate:
    :param tree_alpha:
    :param tree_save_path:
    :param determine_orientation:
    :param randomize_orientation:
    :param orient_decision_boundary:
    :param load_sim_param_from_pkl:
    :return:
    """
    if not os.path.exists(os.path.split(logfile_path)[0]):
        os.makedirs(os.path.split(logfile_path)[0])
    logger = misc.create_logger('Simulation -> Reconstruction', logging.INFO, outfile=logfile_path)

    if load_sim_from_pkl:
        logger.info(f'Loading simulations from {sim_save_path}...')
        with open(os.path.join(sim_save_path, 'dict_crispr_groups.pkl'), 'rb') as f:
            dict_crispr_groups = pickle.load(f)
        # with open(os.path.join(sim_save_path, 'dict_simulation_trees.pkl'), 'rb') as f:
        #     ls_sim_m = pickle.load(f)
        df_sim_protocol = pd.read_pickle(os.path.join(sim_save_path, '_'.join(['0_sim_protocol']) + '.pkl'))
    else:
        logger.info('Starting simulations...')
        df_sim_protocol, dict_crispr_groups, ls_sim_m = run_simulation(sim_parameter_dict, save_path=sim_save_path,
                                                                       plot_tree=sim_plot_tree,
                                                                       logfile_path=sim_logfile_path,
                                                                       sim_as_rec=sim_as_rec,
                                                                       randomize_orientation=randomize_orientation,
                                                                       load_sim_param_from_pkl=load_sim_param_from_pkl,
                                                                       )
        # Alignment
        if sim_as_rec or sim_simple_alignment:
            logger.info('Simulation is used as Reconstruction or sim_simple_alignment was chosen -> '
                        'alignment is done based on simulated top_order.')
            dict_crispr_groups = align_crispr_groups_for_sim_as_rec(dict_crispr_groups)
        else:
            logger.info('Starting alignment...')
            dict_crispr_groups = mafft_integration.align_crispr_groups(WORK_PATH, dict_crispr_groups,
                                                                       mafft_options=alignment_options,
                                                                       logger=logger,)

        if sim_save_path:
            with open(os.path.join(sim_save_path, 'dict_crispr_groups.pkl'), 'wb') as f:
                pickle.dump(dict_crispr_groups, f)

    # If we want to use simulated parameters for grouping.
    if group_by is not None:
        if isinstance(group_by, str):
            if group_by.split('|')[0] == 'sim':
                group_by = df_sim_protocol.set_index('name')[''.join(group_by.split('|')[1:])]
    logger.info('Starting reconstruction...')

    if not os.path.exists(os.path.join(rec_save_path, '0_forward')):
        os.makedirs(os.path.join(rec_save_path, '0_forward'))

    dict_crispr_groups_for_reverse = copy.deepcopy(dict_crispr_groups) if determine_orientation else {}

    (df_rec_protocol, df_rec_protocol_boring, (
        dict_data_for_group_by, lh_fct),
     dict_provided_numbering, dict_trees_forward,
     dict_provided_aligned_arrays,
     dict_provided_duplicated_spacers) = run_reconstruction(
        rec_parameter_dict,
        dict_crispr_groups,
        save_path=os.path.join(
            rec_save_path,
            '0_forward'),
        plot_tree=plot_tree,
        logfile_path=logfile_path,
        do_show=do_show,
        combine_non_unique_arrays=False,
        tree_path=os.path.join(
            sim_save_path,
            'dict_nwk_trees.pkl'),
        use_provided_trees=use_sim_tree,
        tree_save_path=os.path.join(
            tree_save_path,
            '0_forward') if tree_save_path is not None else None,
        sim_as_rec=sim_as_rec,
        hide_unobserved_spacers=hide_unobserved_spacers,
        selection_fct=selection_fct,
        plot_order=plot_order,
        significance_level=significance_level,
        finite_bl_at_root=finite_bl_at_root,
        logger=logger,
        group_by=group_by,
        extend_branches=extend_branches,
        tree_distance_function=tree_distance_function,
        tree_gain_rate=tree_gain_rate,
        tree_loss_rate=tree_loss_rate,
        tree_alpha=tree_alpha,
        alpha_bias_correction=alpha_bias_correction,
        rho_bias_correction=rho_bias_correction,
    )

    dict_trees = dict_trees_forward
    df_rec_protocol = df_rec_protocol.set_index('name')
    if determine_orientation:
        if not os.path.exists(os.path.join(rec_save_path, '0_reversed')):
            os.makedirs(os.path.join(rec_save_path, '0_reversed'))

        for gn, crispr_group in dict_crispr_groups_for_reverse.items():
            for n, ca in crispr_group.crispr_dict.items():
                ca.spacer_array = ca.spacer_array[::-1]
        dict_provided_aligned_arrays = {k: v[v.columns[::-1]] for k, v in dict_provided_aligned_arrays.items()}

        (df_rec_protocol_reversed,
         df_rec_protocol_reversed_boring,
         (dict_data_for_group_by_reversed, lh_fct), _, dict_tree_reversed,
         _, _) = run_reconstruction(rec_parameter_dict,
                                    dict_crispr_groups_for_reverse,
                                    save_path=os.path.join(
                                        rec_save_path,
                                        '0_reversed'),
                                    plot_tree=plot_tree,
                                    logfile_path=logfile_path,
                                    do_show=do_show,
                                    combine_non_unique_arrays=False,
                                    tree_path=os.path.join(
                                        sim_save_path,
                                        'dict_nwk_trees.pkl'),
                                    use_provided_trees=use_sim_tree,
                                    tree_save_path=os.path.join(
                                        tree_save_path,
                                        '0_reversed') if tree_save_path is not None else None,
                                    sim_as_rec=sim_as_rec,
                                    hide_unobserved_spacers=hide_unobserved_spacers,
                                    selection_fct=selection_fct,
                                    plot_order=plot_order,
                                    significance_level=significance_level,
                                    finite_bl_at_root=finite_bl_at_root,
                                    logger=logger,
                                    group_by=group_by,
                                    extend_branches=extend_branches,
                                    tree_distance_function=tree_distance_function,
                                    tree_gain_rate=tree_gain_rate,
                                    tree_loss_rate=tree_loss_rate,
                                    tree_alpha=tree_alpha,
                                    alpha_bias_correction=alpha_bias_correction,
                                    rho_bias_correction=rho_bias_correction,
                                    provided_aligned_arrays=dict_provided_aligned_arrays,
                                    provided_dict_duplicated_spacers=dict_provided_duplicated_spacers,
                                    )
        dict_trees = {'forward': dict_trees_forward, 'reversed': dict_tree_reversed}
        df_rec_protocol_reversed = df_rec_protocol_reversed.set_index('name')
        df_rec_protocol_reversed_boring = df_rec_protocol_reversed_boring.set_index('name')
        df_rec_protocol, df_oriented_rec_protocol, dict_trees = orientation_tools.compare_likelihoods_for_orientation(
            df_rec_protocol,
            df_rec_protocol_reversed,
            df_rec_protocol_boring,
            df_rec_protocol_reversed_boring,
            decision_boundary=orient_decision_boundary, dict_trees=dict_trees)

        if group_by is not None:
            logger.info(f'Grouping trees for parameter estimation by: {group_by}')
            if isinstance(group_by, int):
                val_group_by = []
                for i in range(df_oriented_rec_protocol.shape[0] // group_by):
                    val_group_by += [i] * group_by
                val_group_by += [max(val_group_by)] * (df_oriented_rec_protocol.shape[0] % group_by)
            elif isinstance(group_by, pd.Series):
                val_group_by = [group_by[idx] for idx in df_oriented_rec_protocol.index]
            else:
                val_group_by = df_oriented_rec_protocol[group_by]
            ls_data_for_lh_ratio = []
            filtered_val_group_by = []
            for i, idx in enumerate(df_oriented_rec_protocol.index):
                if df_oriented_rec_protocol.loc[idx, 'our predicted orientation'] == 'Reverse':
                    if idx in dict_data_for_group_by_reversed:
                        ls_data_for_lh_ratio.append(dict_data_for_group_by_reversed[idx])
                        filtered_val_group_by.append(val_group_by[i])
                else:
                    if idx in dict_data_for_group_by:
                        ls_data_for_lh_ratio.append(dict_data_for_group_by[idx])
                        filtered_val_group_by.append(val_group_by[i])

                # ls_data_for_lh_ratio.append(dict_data_for_group_by_reversed[idx] \
                #                                 if df_oriented_rec_protocol.loc[
                #                                        idx, 'our predicted orientation'] == 'Reverse' \
                #                                 else dict_data_for_group_by[idx])
            # logger.info(f'ls_data_for_lh_ratio: {ls_data_for_lh_ratio}')
            # logger.info(f'filtered_val_group_by: {filtered_val_group_by}')
            group_names, groups_lh_0, groups_lh_1, \
                groups_result_0, groups_result_1 = pooling_for_parameter_estimation(ls_data_for_lh_ratio,
                                                                                    group_by=filtered_val_group_by,
                                                                                    give_lh_fct=lh_fct)

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
            df_group_protocol.to_csv(os.path.join(save_path, '0_final_group_by_protocol.csv'))

        df_rec_protocol.to_csv(os.path.join(rec_save_path, '_'.join(['0_protocol_orientation']) + '.csv'))
        df_rec_protocol.to_pickle(os.path.join(rec_save_path, '_'.join(['0_protocol_orientation']) + '.pkl'))
        df_oriented_rec_protocol.to_csv(os.path.join(rec_save_path, '0_final_oriented_protocol.csv'))
        df_oriented_rec_protocol.to_pickle(os.path.join(rec_save_path, '0_final_oriented_protocol.pkl'))

    # trying concat for now, might want to merge/join, especially, if keys overlap
    df_sim_protocol = df_sim_protocol.set_index('name')
    if not determine_orientation:
        if not df_rec_protocol_boring.empty:
            df_rec_protocol_boring = df_rec_protocol_boring.set_index('name')
        df_rec_protocol = pd.concat([df_rec_protocol, df_rec_protocol_boring], axis=0)

    df_rec_protocol = df_oriented_rec_protocol if determine_orientation else df_rec_protocol
    df_sim_rec_protocol = pd.concat([df_sim_protocol, df_rec_protocol], axis=1, join='inner')

    protocol_path = os.path.join(save_path, '_'.join(['0_sim_rec_protocol']))
    logger.info(f'Finished simulation and reconstruction, saving to {protocol_path}')
    df_sim_rec_protocol.to_csv(os.path.join(save_path, '_'.join(['0_sim_rec_protocol']) + '.csv'))
    df_sim_rec_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_sim_rec_protocol']) + '.pkl'))

    if tree_save_path is not None:
        with open(os.path.join(tree_save_path, 'final_dict_nwk_trees.pkl'), 'wb') as f:
            pickle.dump(dict_trees, f)

    # if plot_tree:
    #     for crispr_group in dict_crispr_groups.values():
    #         ls_png = [os.path.join(sim_save_path, crispr_group.name + '_sim.pdf'),
    #                   os.path.join(rec_save_path, crispr_group.name + '_full_rec.pdf')]
    # misc.merge_pdf(ls_png, os.path.join(save_path, crispr_group.name + '_all.pdf'))


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
                       tree_distance_function='breakpoint',
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
                       ):
    """

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
    # rec_m_save_file = open(os.path.join(save_path, '0_rec_m_saved.pkl'), 'wb')
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
        # for key, val in crispr_group.crispr_dict.items():
        #     print(key)
        #     print(val.spacer_array)

        skip, reason = multiple_selection_fcts(crispr_group,
                                               selection_criteria=selection_fct[0] if selection_fct else None,
                                               selection_parameters=selection_fct[1] if selection_fct else None)

        if skip:
            logger.info('Skipped\n for reason: %s', reason)
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

            logger.info(f'Arrays that were combined to one array {dict_combined_arrays}')
            logger.info(f'Number of arrays before: {prev_len} /  after: {len(ls_arrays)}')
            # Be aware: Tree generation or alignment produces an error for single arrays (which are useless anyway)
            if len(ls_arrays) <= 2:
                logger.info('Skipped because arrays were combined. Too few arrays after combining uniques.\n')
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
                logger.info('Skipped\n for reason: %s', reason)
                ls_skipped_protocol.append({'name': crispr_group.name, 'repeat': crispr_group.repeat,
                                            'reason': reason})
                continue
            tree = import_data.load_single_tree_from_string(dict_trees[crispr_group.name] + '\n')
        else:
            # old... If I want support for single file nwk trees, I should do it differently.
            # if tree_path is not None:
            #     logger.info(f'Loading single tree from {tree_path} ...')
            #     tree = import_data.load_single_tree(os.path.join(tree_path, crispr_group.name + '_tree.nwk'))
            # else:
            if lh_fct is None or lh_fct.lower() == 'simple':
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
                        f'provided_lh_fct: {lh_fct} ...')
            tree = construct_tree(ls_array_names, ls_arrays, crispr_group.name,
                                  logger=logger, distance_fct=tree_distance_function,
                                  gain_rate=tree_gain_rate, loss_rate=tree_loss_rate, alpha=tree_alpha,
                                  provided_lh_fct=give_lh_fct,
                                  tree_save_path=None)
        tree = AdvancedTree(tree, True, model_name=crispr_group.name)
        new_dict_trees[crispr_group.name] = tree.format('newick')
        # print(ls_array_names)
        if core_genome_trees:
            # Bio.Phylo.draw_ascii(tree)
            tree, crispr_group = tree_handling(tree, crispr_group, name_inner_nodes=True, bl_eps=0)
            ls_arrays = [crispr_array.spacer_array for crispr_array in crispr_group.crispr_dict.values()]
            ls_array_names = list(crispr_group.crispr_dict.keys())
            if len(ls_arrays) <= 2:
                logger.info('Skipped because too few arrays after pruning.\n')
                ls_skipped_protocol.append({'name': crispr_group.name, 'repeat': crispr_group.repeat,
                                            'reason': 'Skipped because too few arrays after pruning.'})
                continue

        if len(ls_arrays) <= 1:
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
                                   sim_gain_loss_dicts=gain_loss_dicts)

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
                                                  )
            # print('dict_bg_colors', dict_bg_colors)
            # print('Provided numbering: %s' % provided_numbering)
            # print('spacerlabelsnum', spacer_labels_num)

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
        # count_duplications = sum([len(value) for value in rec_duplications.values()])
        # count_rearrangements = sum([len(value) for value in rec_rearrangements.values()])
        # count_reacquisitions = sum([len(value) for value in rec_reacquisitions.values()])
        # count_ind_dups_default = sum([len(value) for value in rec_default_duplications.values()])

        logger.info(f'Length of topological order: {len(rec_m.top_order)}')
        nb_gained_spacers = len(rec_m.top_order)

        nb_rec_insertions = sum([len(val) for val in rec_m.rec_gain_dict.values()])
        nb_rec_deletions = sum([len(val) for val in rec_m.rec_loss_dict.values()])
        # unique_spacers = set()
        # for leaf in rec_m.get_terminals():
        #     unique_spacers = unique_spacers.union(set(leaf.spacers))
        nb_unique_spacers = len(rec_m.top_order) - count_duplications_rearrangement_candidates
        # logger.info(f'All spacer orders: {rec_m.spacer_orders}')
        dict_matrices = rec_m.generate_possible_losses_sets()
        # logger.info(f'Possible multiple losses sets: {dict_matrices}')
        dict_max_length_losses, ls_all_max_lengths, ls_all_max_lengths_norm = rec_m.get_max_length_loss_sets()
        # logger.info(f'Max length multiple losses sets (set, length): {dict_max_length_losses}')
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
        # ln_lh_ratio = 2 * (opt_result_0.fun - opt_result_1.fun)
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
                         'run time': run_time,
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
        # pickle.dump((i, rec_m.get_data_for_lh_ratio(filter_by=False)), rec_m_save_file)
        logger.info(f'Run time of group "{crispr_group.name}": {run_time} s\n')
        # tree_memory_size = sys.getsizeof(rec_m)
        # memory_size_of_function = sys.getsizeof(give_lh_fct)
        # logger.info(f'Memory size of tree: {tree_memory_size}, Memory size lh function: {memory_size_of_function}')
        # local_vars = list(locals().items())
        # for var, obj in local_vars:
        #     print(var, sys.getsizeof(obj))
        # for name, size in sorted(((name, misc.getsize(value)) for name, value in locals().items()),
        #                          key=lambda x: -x[1])[:20]:
        #     print("{:>30}: {:>8}".format(name, misc.sizeof_fmt(size)))
    # rec_m_save_file.close()
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
        # ln_lh_ratio = 2 * (opt_result_0.fun - opt_result_1.fun)
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


def run_simulation(sim_parameter_dict, save_path=None, plot_tree=True, logfile_path=None,
                   do_show=False,
                   sim_as_rec=False,
                   randomize_orientation=False,
                   load_sim_param_from_pkl=None,
                   ):
    """
    :param load_sim_param_from_pkl:
    :param randomize_orientation:
    :param sim_as_rec:
    :param sim_parameter_dict:
    :param save_path:
    :param plot_tree:
    :param logfile_path:
    :param do_show:
    :return:
    """
    if not os.path.exists(os.path.split(logfile_path)[0]):
        os.makedirs(os.path.split(logfile_path)[0])
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    logger = misc.create_logger('simulation', logging.INFO, outfile=logfile_path)
    logger.info('Running Simulation of parameter dictionary: %s', sim_parameter_dict)
    run_timer = misc.RunTimer()

    if sim_parameter_dict.get('import_path', False):
        ls_tree_names, ls_trees, ls_gtr_gain_rate, ls_gtr_loss_rate, ls_deletion_model, ls_deletion_model_parameter, \
            ls_sim_seed, ls_repeat_runs = expand_sim_parameters_based_on_import(sim_parameter_dict)
    else:
        ls_tree_names, ls_trees, ls_gtr_gain_rate, ls_gtr_loss_rate, ls_deletion_model, ls_deletion_model_parameter, \
            ls_sim_seed = expand_sim_parameters(sim_parameter_dict)
        ls_repeat_runs = [0] * len(ls_tree_names)

    ls_dict_protocol = []
    dict_sim_m = {}
    dict_crispr_groups = {}
    dict_trees = {}
    ls_skipped_protocol = []

    for i, (tree, name) in enumerate(zip(ls_trees, ls_tree_names)):
        repeat_run = ls_repeat_runs[i]
        deletion_model = ls_deletion_model[i]
        gtr_gain_rate = ls_gtr_gain_rate[i]
        gtr_loss_rate = ls_gtr_loss_rate[i]
        sim_seed = ls_sim_seed[i]

        logger.info('Working on ' + name)
        loss_rate = ls_gtr_loss_rate[i]
        gain_rate = ls_gtr_gain_rate[i]

        deletion_model = ls_deletion_model[i] if isinstance(deletion_model, list) else deletion_model
        deletion_model_parameter = ls_deletion_model_parameter[i]
        logger.info(f'Gain rate: {gain_rate}, Loss rate: {loss_rate}')
        logger.info(f'Deletion model: {deletion_model}')

        if gain_rate == 0 or loss_rate == 0:
            logger.info('Skipped\n')
            ls_skipped_protocol.append({'name': name, 'reason': 'gain or loss rate is 0'})
            continue

        # give names to nodes, if necessary
        n_i = 0
        for node in tree.get_nonterminals():
            if node.name is None:
                node.name = 'in_' + str(n_i)
                n_i += 1

        m = SimulationTree(deletion_model, sim_parameter_dict['event_sim'], gtr_gain_rate, gtr_loss_rate,
                           deletion_model_parameter,
                           os.path.join(save_path), tree, True, model_name=name,
                           blm_rate_correction_last_spacer=sim_parameter_dict.get('blm_rate_correction_last_spacer',
                                                                                  True),
                           blm_two_directional=sim_parameter_dict.get('blm_two_directional', False), )

        tree_height = m.distance_to_leafs()
        tree_length = m.total_branch_length()
        tree_leafs = m.nb_leafs()
        logger.info('Tree parameters: Distance to leaf: %s ; tree length: %s ; nb of leafs: %s' % (tree_height,
                                                                                                   tree_length,
                                                                                                   tree_leafs))

        m.simulation_on_tree(start_spacers='sample', seed=sim_seed)

        gained_spacers = m.next_spacer - 1
        avg_array_len = m.avg_length_leaf_spacers()
        logger.info(f'Number of gained spacers: {gained_spacers}')
        logger.info(f'Average length of leaf spacers: {avg_array_len}')
        never_seen_spacers, nb_never_seen_spacers = m.get_never_seen_spacers()
        logger.info(f'Nb of not observed spacers: {len(never_seen_spacers)}')
        if sim_parameter_dict['event_sim']:
            ls_evt_times, ls_evt_losses, ls_evt_nb_ex_spacers, \
                ls_evt_loss_pos, ls_evt_rel_loss_pos = m.get_sim_rel_loss_pos()
        else:
            ls_evt_rel_loss_pos, ls_evt_nb_ex_spacers = [], []

        (ls_branch_based_rel_loss_pos, ls_branch_based_losses,
         ls_branch_based_nb_existent_spacers) = m.get_relative_loss_positions_branch_based()

        if plot_tree:
            m.visualize_tree(name=name, do_show=do_show, path=save_path, )

        dict_sim_m[name] = m
        dict_crispr_groups[name] = m.get_data_as_crispr_group(name=name, only_terminals=not sim_as_rec)
        dict_trees[name] = AdvancedTree(tree, True).format(fmt='newick')

        (fes_m_presence, les_m_presence, l_s_m_presence, f_s_m_presence), \
            (global_fes_m_presence, global_les_m_presence, global_l_s_presence) = m.get_relative_loss_stats()

        sim_orientation = 1
        if randomize_orientation:
            if random.randint(1, 2) == 2:
                sim_orientation = 2
                for n, ca in dict_crispr_groups[name].crispr_dict.items():
                    ca.spacer_array = ca.spacer_array[::-1]
                    ca.orientation = sim_orientation

        logger.info(f'Run time of experiment "{name}": {run_timer.time_from_last_checkpoint()}\n')

        dict_protocol = {'name': name,
                         'sim. gain rate': gain_rate, 'sim. loss rate': loss_rate,
                         'sim loss model': deletion_model,
                         'sim avg. array length': avg_array_len,
                         'sim leaf array lengths': m.get_array_lengths(),
                         'loss model param.': deletion_model_parameter,
                         'not observed spacers': never_seen_spacers,
                         'nb not observed spacers': nb_never_seen_spacers,
                         'sim. relative loss positions': ls_evt_rel_loss_pos,
                         'sim. branch based relative loss positions': ls_branch_based_rel_loss_pos,
                         'sim branch based nb of existent spacers': ls_branch_based_nb_existent_spacers,
                         'sim nb of existent spacers': ls_evt_nb_ex_spacers,
                         'sim seed': sim_seed,
                         'tree seed of each tree batch': sim_parameter_dict.get('tree_seed', None),
                         'sim orientation': sim_orientation,
                         'repeat run': repeat_run,
                         'sim fes +m presence': fes_m_presence, 'sim les -m presence': les_m_presence,
                         'sim ls -m presence': l_s_m_presence,
                         'sim fs +m presence': f_s_m_presence,
                         'sim global fes +m presence': global_fes_m_presence,
                         'sim global les -m presence': global_les_m_presence,
                         'sim global ls -m presence': global_l_s_presence,
                         }

        ls_dict_protocol.append(dict_protocol)

    logger.info(f'Runtime of all experiments: {run_timer.time_from_start()}')

    if save_path:
        logger.info(f'Saving model to {save_path}')
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        # with open(os.path.join(save_path, 'dict_simulation_trees.pkl'), 'wb') as f:
        #     pickle.dump(dict_sim_m, f)
        with open(os.path.join(save_path, 'dict_nwk_trees.pkl'), 'wb') as f:
            pickle.dump(dict_trees, f)

    df_protocol = pd.DataFrame(ls_dict_protocol)
    df_protocol.to_csv(os.path.join(save_path, '_'.join(['0_sim_protocol']) + '.csv'))
    df_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_sim_protocol']) + '.pkl'))
    df_skipped = pd.DataFrame(ls_skipped_protocol)
    df_skipped.to_csv(os.path.join(save_path, '_'.join(['0_skipped_protocol']) + '.csv'))

    return df_protocol, dict_crispr_groups, dict_sim_m
