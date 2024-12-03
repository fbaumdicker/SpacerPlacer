import logging
import os
import pandas as pd
import pickle
import random

from sp_model.helpers import misc
from sp_model.experiments_tools import expand_sim_parameters, expand_sim_parameters_based_on_import
from sp_model.simulation_tree import SimulationTree
from sp_model.data_classes.advanced_tree import AdvancedTree


def run_simulation(sim_parameter_dict, save_path=None, plot_tree=True, logfile_path=None,
                   do_show=False,
                   sim_as_rec=False,
                   randomize_orientation=False,
                   ):
    """
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
                         'repeat main': repeat_run,
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
        with open(os.path.join(save_path, 'dict_nwk_trees.pkl'), 'wb') as f:
            pickle.dump(dict_trees, f)

    df_protocol = pd.DataFrame(ls_dict_protocol)
    df_protocol.to_csv(os.path.join(save_path, '_'.join(['0_sim_protocol']) + '.csv'))
    df_protocol.to_pickle(os.path.join(save_path, '_'.join(['0_sim_protocol']) + '.pkl'))
    df_skipped = pd.DataFrame(ls_skipped_protocol)
    df_skipped.to_csv(os.path.join(save_path, '_'.join(['0_skipped_protocol']) + '.csv'))

    return df_protocol, dict_crispr_groups, dict_sim_m
