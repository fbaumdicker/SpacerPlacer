import copy

import pandas as pd
import numpy as np
import os
import pickle

from sp_model.helpers import misc, import_data
from sp_model import model_tools
from sp_model.coalescent_tree_sim import mult_rnd_coalescent_tree, rnd_coalescent_tree
from sp_model import distance_tree_constructor
from sp_model import mafft_integration
from sp_model.data_classes.crisprdata import CRISPRGroup, CRISPRArray, Species
from sp_model.data_classes.advanced_tree import AdvancedTree

GROUPS_TO_EXCLUDE = ['g_39', 'g_284', 'g_690', 'g_140', 'g_719', 'g_931', 'g_598',
                     'g_942']


def load_align_pickled_data(path, mafft_options=None, exclude_files=None, save_path=None):
    """
    :param path:
    :param mafft_options:
    :param exclude_files:
    :param save_path:
    :return:
    """
    ls_files = []
    for entry in sorted(os.listdir(path)):
        if os.path.isdir(os.path.join(path, entry)):
            continue
        ls_files.append(entry)
    dict_crispr_groups = {}
    if exclude_files:
        ls_files = [file for file in ls_files if file not in exclude_files]
    for file in ls_files:
        ls_arrays, ls_names, head = mafft_integration.pickled_group_to_mafft(path, [file],
                                                                             options=mafft_options)
        ls_crispr_arrays = []
        for name, array in zip(ls_names, ls_arrays):
            crispr_array = CRISPRArray(name, '', None, 4, head[1],
                                       'chromosome', 'Bacteria', head[0], array, head[2],
                                       cas_genes=[], all_cas_types=[], all_cas_genes=[],
                                       species=Species('', 'Bacteria', 'sim', 'sim', 'sim', 'sim',
                                                                      'sim', 'sim'),
                                       species_fc='sim')
            ls_crispr_arrays.append(crispr_array)
        crispr_group = CRISPRGroup(head[1], ls_crispr_arrays, name=file)
        crispr_group.convert_arrays_from_mafft_fmt()
        dict_crispr_groups[file] = crispr_group

    if save_path:
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        with open(os.path.join(save_path, '0_aligned_crisprgroups.pkl'), 'wb') as f:
            pickle.dump(dict_crispr_groups, f)
    return dict_crispr_groups


def give_count_duplicates_in_list(ls):
    """
    Returns a list of all duplicates in the list
    :param ls:
    :return:
    """
    dict_duplicates = {}
    ls_unique = []
    for item in ls:
        if item in ls_unique:
            dict_duplicates[item] = dict_duplicates.get(item, 1) + 1
        else:
            ls_unique.append(item)
    return dict_duplicates


def name_inner_nodes_if_unnamed(tree):
    n_i = 0
    for node in tree.get_nonterminals():
        if node.name is None:
            node.name = 'in_' + str(n_i)
            n_i += 1
    return tree


def tree_handling(tree, crispr_group, name_inner_nodes=False, bl_eps=0):
    dict_duplicates = give_count_duplicates_in_list(crispr_group.get_ls_acc_nums())
    duplicate_to_ls_arrays = {key: crispr_group.get_ls_arrays_by_acc_num(key) for key in dict_duplicates.keys()}
    # acc_num: name of arrays that need to be pruned
    dict_to_prune = {}
    # acc_num to names of arrays that need to use the same leaf... (or combine them?)
    dict_acc_num_ls_duplicates = {}
    for key, ls_ca in duplicate_to_ls_arrays.items():
        ls_ca_w_o_gaps = [[a for a in x.spacer_array if a != '-'] for x in ls_ca]
        if all([ls_ca_w_o_gaps[0] == x for x in ls_ca_w_o_gaps]):
            dict_acc_num_ls_duplicates[key] = [x.name for x in ls_ca]
        else:
            dict_to_prune[key] = [x.name for x in ls_ca]

    ls_crispr_arrays_to_delete = []
    for name, ca in crispr_group.crispr_dict.items():
        for leaf in tree.get_terminals():
            if leaf.name == ca.acc_num:
                if leaf.name in dict_to_prune.keys():
                    # print('pruning', leaf.name)
                    if len(list(tree.get_terminals())) > 1:
                        tree.prune(leaf)
                    ls_crispr_arrays_to_delete.append(ca.acc_num)
                # This splits the tree leafs for equal CRISPR arrays in the same strain...
                elif leaf.name in dict_acc_num_ls_duplicates.keys():
                    leaf.split(n=len(dict_acc_num_ls_duplicates[leaf.name]), branch_length=bl_eps)
                    for child in leaf.clades:
                        child.name = dict_acc_num_ls_duplicates[leaf.name].pop()
                else:
                    leaf.name = name
    for acc_num in ls_crispr_arrays_to_delete:
        crispr_group.drop_array_by_acc_num(acc_num)
    # Naming inner nodes...
    if name_inner_nodes:
        tree = name_inner_nodes_if_unnamed(tree)
    tree = AdvancedTree(tree, rooted=True, model_name=tree.model_name)
    return tree, crispr_group


def load_align_single_fasta(data_path, work_path, mafft_options=None, group_name='g', logger=None, save_path=None,
                            seed=None):
    dict_fasta = read_fasta(data_path)
    # group metadata format?
    # repeat = 'ABC'

    ls_array_names, ls_spacer_arrays = [], []
    # ls_array_metadata = []
    for key, value in dict_fasta.items():
        split_line = key.split(',')
        ls_array_names.append(split_line[0])
        # ls_array_metadata.append(split_line[1:])
        ls_spacer_arrays.append(value.replace(' ', '').split(','))
    ls_crispr_arrays = []
    for name, array in zip(ls_array_names, ls_spacer_arrays):
        crispr_array = CRISPRArray(name, None, None, None,
                                   None, None, None, None, array,
                                   None,
                                   cas_genes=None, all_cas_types=None, all_cas_genes=None,
                                   species=None,
                                   species_fc=None)
        ls_crispr_arrays.append(crispr_array)
    crispr_group = CRISPRGroup(None, ls_crispr_arrays, name=group_name)
    dict_crispr_group = {crispr_group.name: crispr_group}
    dict_crispr_group = mafft_integration.align_crispr_groups(work_path, dict_crispr_group, mafft_options=mafft_options,
                                                              logger=logger, seed=seed)
    if save_path:
        with open(save_path, 'wb') as f:
            pickle.dump(dict_crispr_group, f)
    return list(dict_crispr_group.values())[0]


def read_fasta(data_path):
    with open(data_path, 'r') as f:
        lines = f.readlines()
        headers = []
        sequences = []
        for line in lines:
            if line[0] == '>':
                headers.append(line[1:].strip())
            else:
                sequences.append(line.strip())
        return {header: sequence for header, sequence in zip(headers, sequences)}


def align_crispr_groups_for_sim_as_rec(dict_crispr_groups):
    for key, crispr_group in dict_crispr_groups.items():
        dict_arrays_as_list = crispr_group.get_arrays_as_lists()
        top_order = crispr_group.top_order
        ls_names = list(dict_arrays_as_list.keys())
        ls_arrays = list(dict_arrays_as_list.values())

        ls_aligned_arrays = []
        for array in ls_arrays:
            aligned_array = [str(s) if str(s) in array else '-' for s in top_order]
            ls_aligned_arrays.append(aligned_array)
        crispr_group.update_spacer_arrays_by_ls_arrays_as_list(ls_names, ls_aligned_arrays)
    return dict_crispr_groups


def filter_fct(dict_protocol):
    """
    :param dict_protocol:
    :return:
    """
    # return True
    return True if dict_protocol['all max block deletion lengths'] else False


def construct_tree(ls_array_names, ls_arrays, group_name, logger=None, distance_fct='breakpoint', tree_save_path=None,
                   gain_rate=None, loss_rate=None, alpha=None, provided_lh_fct=None, tree_construction_method='upgma'):
    """
    Constructs an upgma tree based on different distance functions.
    :param tree_construction_method:
    :param ls_array_names:
    :param ls_arrays:
    :param group_name:
    :param logger:
    :param distance_fct:
    :param tree_save_path:
    :param gain_rate:
    :param loss_rate:
    :param alpha:
    :param provided_lh_fct:
    :return:
    """
    if distance_fct == 'likelihood':
        distance = distance_tree_constructor.LikelihoodDistance(gain_rate, loss_rate, alpha,
                                                                provided_lh_fct=provided_lh_fct, skip_letters={'-'})
    else:
        raise NotImplementedError('Distance function not implemented!')

    tree_constructor = distance_tree_constructor.FixedDistanceTreeConstructor(distance_calculator=distance,
                                                                              method=tree_construction_method)
    # tree_constructor = DistanceTreeConstructor(distance_calculator=distance, method=tree_construction_method)

    if logger is not None:
        logger.info(f'Constructing tree of {group_name}...')
    dm = tree_constructor.distance_calculator.get_distance([ls_array_names, ls_arrays])
    tree = tree_constructor.build_tree([ls_array_names, ls_arrays])


    # for clade in tree.find_clades():
    #     print(clade.name)
    #     print(clade.branch_length)
    if tree_save_path:
        if not os.path.exists(tree_save_path):
            os.makedirs(tree_save_path)
        save_name = os.path.join(tree_save_path, group_name + '_tree.nwk')
        logger.info(f'Tree saved to: {save_name}')
        misc.write_tree_to_file(tree, save_name)
    return tree


def multiple_selection_fcts(array, selection_criteria, selection_parameters=None):
    skip = False
    reason = None
    if not selection_criteria:
        skip = False
    else:
        for selection_criterium, selection_parameter in zip(selection_criteria, selection_parameters):
            skip, reason = selection_fct(array, selection_criterium, selection_parameter=selection_parameter)
            if skip:
                return skip, reason
    return skip, reason


def selection_fct(group, selection_criterium, selection_parameter=None):
    if selection_criterium is None:
        skip = False
    elif selection_criterium == 'maximum_unique_array_number':
        ls_arrays = []
        for array in group.crispr_dict.values():
            if array.spacer_array in ls_arrays:
                continue
            else:
                ls_arrays.append(array.spacer_array)
        skip = True if len(ls_arrays) > selection_parameter else False
    elif selection_criterium == 'maximum_array_number':
        skip = True if len(group.crispr_dict) > selection_parameter else False
    elif selection_criterium == 'minimum_array_number':
        skip = True if len(group.crispr_dict) < selection_parameter else False
    elif selection_criterium == 'empty_arrays':
        skip = True
        for key, array in group.crispr_dict.items():
            if set(array.spacer_array) not in [{'-'}, set()]:
                skip = False
                break
    elif selection_criterium == 'crispr_type':
        if group.cas_type.lower() in selection_parameter:
            skip = False
        else:
            skip = True
    elif selection_criterium == 'selection of groups':
        if group.name in selection_parameter:
            skip = False
        else:
            skip = True
    elif selection_criterium == 'selection not in groups':
        if group.name in selection_parameter:
            skip = True
        else:
            skip = False
    else:
        raise NotImplementedError('Selection function not implemented!')
    reason = selection_criterium
    return skip, reason


def load_arrays_from_pkl(data_path):
    dict_crispr_arrays = import_data.load_crispr_arrays_from_pickle(os.path.join(data_path))
    return dict_crispr_arrays


def expand_sim_parameters_based_on_import(parameter_dict):
    """
    needs import_path for parameters, trees_path for trees, and repeat_runs for number of repeated runs.
    :param parameter_dict:
    :return:
    """
    import_path = parameter_dict['import_path']
    with open(os.path.join(import_path), 'rb') as f:
        df_imported_parameters = pd.read_pickle(f)
        print('imported parameter df head: ', df_imported_parameters.head())
        ls_names = df_imported_parameters.index.values
        ls_names = [name for name in ls_names if name not in GROUPS_TO_EXCLUDE]
        df_imported_parameters = df_imported_parameters.loc[ls_names]
        # df_imported_parameters = df_imported_parameters.set_index('name')

    ls_avg_array_length = list(df_imported_parameters['avg array length'].values)
    if parameter_dict['deletion_model'] == 'independent':
        ls_gtr_loss_rate = list(df_imported_parameters['loss_rate_0'].values)
        ls_deletion_model = ['independent'] * len(ls_gtr_loss_rate)
        ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_length,
                                                                  ls_gtr_loss_rate,
                                                                  'independent'
                                                                  )
        ls_deletion_model_parameter = [None] * len(ls_gtr_loss_rate)
    elif parameter_dict['deletion_model'] == 'block':
        ls_gtr_loss_rate = list(df_imported_parameters['loss_rate_1'].values)
        ls_deletion_model = ['block'] * len(ls_gtr_loss_rate)
        ls_deletion_model_parameter = list(df_imported_parameters['alpha_1'].values)
        ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_length,
                                                                  ls_gtr_loss_rate,
                                                                  'block',
                                                                  ls_deletion_model_parameter,
                                                                  )
    extended_ls_tree_names = []
    extended_ls_trees = []
    extended_ls_gtr_gain_rate = []
    extended_ls_gtr_loss_rate = []
    extended_ls_deletion_model = []
    extended_ls_deletion_model_parameter = []
    extended_ls_repeat_runs = []
    repeat_runs = parameter_dict['repeat_runs']
    extended_ls_sim_seed = list(range(parameter_dict['sim_seed'],
                                      parameter_dict['sim_seed'] + len(ls_names) * repeat_runs))

    nb_leafs_col = df_imported_parameters['nb leafs (after combining)']
    if parameter_dict.get('trees_path', False):
        trees_path = parameter_dict['trees_path']
        with open(trees_path, 'rb') as f:
            dict_trees = pickle.load(f)
    elif parameter_dict.get('tree_sim_per_repeat_run', False):
        ls_tree_sim_seed = [parameter_dict['tree_seed'] + i for i in range(len(ls_names))] \
            if 'tree_seed' in parameter_dict else [None] * len(ls_names)
        dict_trees = {ls_names[i]: rnd_coalescent_tree(nb_leafs_col.iloc[i], rnd_bl=True, seed=ls_tree_sim_seed[i],
                                                       tree_name=ls_names[i]) for i in range(len(ls_names))}
    elif parameter_dict['tree_sim_per_repeat_run']:
        dict_trees = {}
        ls_tree_sim_seed = [parameter_dict['tree_seed'] + i for i in range(repeat_runs * len(ls_names))] \
            if 'tree_seed' in parameter_dict else [None] * (repeat_runs * len(ls_names))
    else:
        dict_trees = None

    for r in range(repeat_runs):
        new_ls_tree_names = [name + '-' + str(r) for name in ls_names]
        extended_ls_tree_names += new_ls_tree_names
        if parameter_dict.get('trees_path', False) or parameter_dict.get('tree_sim_per_repeat_run', False):
            new_ls_trees = [copy.copy(dict_trees[name]) for name in ls_names]
            new_ls_trees = [import_data.load_single_tree_from_string(tree) for tree in new_ls_trees]
        elif parameter_dict.get('tree_sim_per_repeat_run', False):
            new_ls_trees = [rnd_coalescent_tree(nb_leafs_col.loc[name],
                                                rnd_bl=True,
                                                seed=ls_tree_sim_seed[r * len(ls_names) + i],
                                                tree_name=new_ls_tree_names[i])
                            for i, name in enumerate(ls_names)]
        for i, t in enumerate(new_ls_trees):
            t.name = new_ls_tree_names[i]
        extended_ls_repeat_runs += [r] * len(new_ls_tree_names)
        extended_ls_trees += new_ls_trees
        extended_ls_gtr_gain_rate += ls_gtr_gain_rate
        extended_ls_gtr_loss_rate += ls_gtr_loss_rate
        extended_ls_deletion_model += ls_deletion_model
        extended_ls_deletion_model_parameter += ls_deletion_model_parameter

    return extended_ls_tree_names, extended_ls_trees, extended_ls_gtr_gain_rate, extended_ls_gtr_loss_rate, \
        extended_ls_deletion_model, extended_ls_deletion_model_parameter, extended_ls_sim_seed, extended_ls_repeat_runs


def expand_sim_parameters(parameter_dict):
    """
    Should include: nb_trees, nb_leafs, deletion_model (string or list), gtr_loss_rate, deletion_model
    tree_seed for consistency
    deletion_model_parameter if used
    parameters_path if used
    lr_dependent_on_length_alpha desired_length if used then no gtr_gain_rate
    trees_path if used then no rnd_bl
    :param parameter_dict:
    :return: dictionary with trees, gtr_gain/loss_rate, deletion_model, deletion_model_parameter
    """
    nb_trees = parameter_dict['nb_trees']
    nb_leafs = parameter_dict['nb_leafs']
    ls_tree_names = None

    sim_seed = parameter_dict.get('sim_seed', None)

    if parameter_dict.get('parameters_path', False):
        parameters_path = parameter_dict['parameters_path']
        df_param = pd.read_pickle(parameters_path)
        ls_tree_names = list(df_param['name'])
        ls_avg_array_len = list(df_param['avg array length'])
        if parameter_dict['deletion_model'] == 'independent':
            ls_gtr_loss_rate = list(df_param['loss_rate_0'])
            ls_deletion_model = ['independent'] * len(ls_gtr_loss_rate)
            ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_len,
                                                                      ls_gtr_loss_rate,
                                                                      'independent'
                                                                      )
            ls_deletion_model_parameter = [None] * len(ls_gtr_gain_rate)
        elif parameter_dict['deletion_model'] == 'block':
            ls_gtr_loss_rate = list(df_param['loss_rate_1'])
            ls_deletion_model_parameter = list(df_param['alpha_1'])
            ls_deletion_model = ['block'] * len(ls_gtr_loss_rate)
            ls_gtr_gain_rate = misc.estimate_gr_based_on_estimated_lr(ls_avg_array_len,
                                                                      ls_gtr_loss_rate,
                                                                      'block',
                                                                      ls_deletion_model_parameter,
                                                                      )
        else:
            raise NotImplementedError('Model not implemented for importing parameters!')
        ls_sim_seed = [None] * len(ls_gtr_loss_rate) if sim_seed is None else list(
            range(sim_seed, sim_seed + len(ls_gtr_loss_rate)))
    else:
        deletion_model = parameter_dict['deletion_model']
        deletion_model_parameter = parameter_dict.get('deletion_model_parameter', None)
        if parameter_dict['lr_based_on_gr_avg_array_length']:
            gtr_gain_rate = parameter_dict['gtr_gain_rate']
            gtr_loss_rate = misc.keep_len_constant(gtr_gain_rate, parameter_dict['avg_array_length'], deletion_model,
                                                   ls_alpha=deletion_model_parameter)
        else:
            gtr_gain_rate = parameter_dict['gtr_gain_rate']
            gtr_loss_rate = parameter_dict['gtr_loss_rate']
        nb_trees = parameter_dict['nb_trees']
        nb_leafs = parameter_dict['nb_leafs']
        ls_gtr_gain_rate = []
        ls_gtr_loss_rate = []
        ls_deletion_model = []
        ls_deletion_model_parameter = []
        ls_sim_seed = []
        for i, (t, l) in enumerate(zip(nb_trees, nb_leafs)):
            if isinstance(gtr_gain_rate, list):
                ls_gtr_gain_rate += list(np.repeat(gtr_gain_rate[i], t))
                ls_gtr_loss_rate += list(np.repeat(gtr_loss_rate[i], t))
            else:
                ls_gtr_gain_rate += list(np.repeat(gtr_gain_rate, t))
                ls_gtr_loss_rate += list(np.repeat(gtr_loss_rate, t))
            ls_deletion_model += list(np.repeat(deletion_model[i], t))
            ls_sim_seed += [None] * t if sim_seed is None else list(range(sim_seed, sim_seed + t))
            if deletion_model_parameter:
                ls_deletion_model_parameter += list(np.repeat(deletion_model_parameter[i], t))
            else:
                ls_deletion_model_parameter += [None] * t

    if parameter_dict.get('trees_path', False):
        if ls_tree_names is None:
            ls_tree_names = parameter_dict('tree_file_names')
        if parameter_dict.get('trees_path', None):
            if parameter_dict.get('trees_format', 'nwk'):
                ls_trees = import_data.load_real_data_tree(os.path.join(parameter_dict['trees_path']),
                                                           [n + '_tree.nwk' for n in ls_tree_names])
            else:
                with open(parameter_dict['trees_path']) as f:
                    ls_trees = pickle.load(f)
        for n, t in zip(ls_tree_names, ls_trees):
            t.name = n
    else:
        ls_trees = []
        ls_tree_names = []
        tree_seed = parameter_dict.get('tree_seed', None)
        for i, (t, l) in enumerate(zip(nb_trees, nb_leafs)):
            ls_trees += mult_rnd_coalescent_tree(t, l, rnd_bl=parameter_dict['rnd_bl'], seed=tree_seed,
                                                 tree_name='t')
        ls_tree_names += ['_'.join([parameter_dict['model_name'], str(i)]) for i, t in enumerate(ls_trees)]

    return ls_tree_names, ls_trees, ls_gtr_gain_rate, ls_gtr_loss_rate, ls_deletion_model, \
        ls_deletion_model_parameter, ls_sim_seed


def pooling_for_parameter_estimation(ls_data, group_by=None, give_lh_fct=None):
    values = np.unique(group_by)
    ls_grouped_idx = []
    for val in values:
        group_idx = []
        for i, x in enumerate(group_by):
            if x == val:
                group_idx.append(i)
        ls_grouped_idx.append(group_idx)
    ls_lh_0, ls_lh_1, ls_opt_result_0, ls_opt_result_1 = [], [], [], []
    for group_idx in ls_grouped_idx:
        group_data = [ls_data[idx] for idx in group_idx]
        lh_0, lh_1, opt_result_0, opt_result_1 = model_tools.compute_lh_ratio_of_multiple_trees(group_data,
                                                                                                method='Nelder-Mead',
                                                                                                filter_by=False,
                                                                                                lambdifyed_lh_fct=give_lh_fct)
        ls_lh_0.append(lh_0)
        ls_lh_1.append(lh_1)
        ls_opt_result_0.append(opt_result_0)
        ls_opt_result_1.append(opt_result_1)
    return values, ls_lh_0, ls_lh_1, ls_opt_result_0, ls_opt_result_1
