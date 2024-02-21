import sys
import os
import pandas as pd

import experiments_run

MAFFT_OPTIONS = ['--text',
                 # '--genafpair',
                 '--localpair',
                 # '--sp',
                 # '--globalpair',
                 # '--hybridpair',
                 # '--leavegappyregion',
                 '--maxiterate', '1000',
                 # '--noscore',
                 '--lep', '0',
                 '--op', '0',
                 '--lop', '0',
                 '--lep', '0',
                 '--lexp', '0',
                 '--LOP', '0',
                 '--LEXP', '0',
                 '--randomseed', '2357',
                 '--quiet',
                 '--thread', '-1',  # multithreading -1 -> automatically choose
                 # '--threadit', '0',  # threaded iterative alignment has some randomness (this stops this)
                 ]

TREE_GENERATION_PARAMETERS = {'old_gain_rate': 0.9724,  # based on omer data (without boring data; old and new lh fct)!
                              'old_loss_rate': 0.1765,
                              'old_alpha': 5.5092,
                              'new_gain_rate': 0.6108,
                              'new_loss_rate': 0.1830,
                              'new_alpha': 3.3377,
                              }


def run_omer_data():
    finite_bl_at_root = True
    extend_branches = 0.01  # 0.00001
    group_by = 'crispr type'
    rec_parameter_dict = {'model': 'gtr', 'gain_rate': 1, 'loss_rate': 0.1,
                          # 'give_lh_fct': os.path.join('data', '0_lh', 'death_lh_up_to_100_lambdified.pkl'),
                          'give_lh_fct': None
                          }

    tree_distance_function = 'breakpoint'
    if rec_parameter_dict['give_lh_fct'] is None:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['old_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['old_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['old_alpha']
    else:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['new_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['new_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['new_alpha']

    # data_path = os.path.join('data', 'omer_data')
    data_path = os.path.join('data', 'omer_aligned_data', '0_aligned_crisprgroups.pkl')
    alignment_save_path = os.path.join('data', 'omer_aligned_data')
    selection_fct = (['crispr_type', 'minimum_array_number', 'selection not in groups'],
                     [['cas-i-e', 'cas-i-f', 'cas-ii-a', 'cas-i-c', 'cas-ii-c', 'cas-iii-a'], 2, ['group_561']],
                     )
    # selection_fct = (['selection of groups'], [['group_561']])
    result_folder = os.path.join('0_result_folder', 'omer_data')
    save_path = os.path.join(result_folder, 'gb_crispr_type_old_lh_fct_fin_bl')
    logfile_path = os.path.join(save_path, '0_logger.log')
    combine_non_unique_arrays = True
    tree_path = None
    tree_save_path = None
    hide_unobserved_spacers = False
    plot_order = True
    plot_tree = True
    significance_value = 0.05

    alignment_options_dict = {'mafft_options': MAFFT_OPTIONS,
                              'exclude_files': ['group_707', 'group_672', 'group_711', 'group_710', 'group_1027',
                                                'group_357',
                                                ],
                              'save_path': alignment_save_path}

    experiments_run.run_pickled_data(rec_parameter_dict, data_path=data_path, save_path=save_path,
                                     plot_tree=plot_tree, logfile_path=logfile_path, do_show=False,
                                     combine_non_unique_arrays=combine_non_unique_arrays, tree_path=tree_path,
                                     tree_save_path=tree_save_path, hide_unobserved_spacers=hide_unobserved_spacers,
                                     selection_fct=selection_fct,
                                     plot_order=plot_order, significance_level=significance_value,
                                     finite_bl_at_root=finite_bl_at_root,
                                     alignment_options_dict=alignment_options_dict,
                                     group_by=group_by,
                                     extend_branches=extend_branches,
                                     tree_distance_function=tree_distance_function,
                                     tree_gain_rate=tree_gain_rate,
                                     tree_loss_rate=tree_loss_rate,
                                     tree_alpha=tree_alpha,
                                     )


def run_crispr_db_data():
    determine_orientation = True
    orient_boundary = 10
    finite_bl_at_root = True
    extend_branches = 0.00001
    group_by = 'crispr type'
    rec_parameter_dict = {'model': 'gtr', 'gain_rate': 1, 'loss_rate': 0.1,
                          # 'give_lh_fct': os.path.join('data', '0_lh', 'death_lh_up_to_100_lambdified.pkl'),
                          'give_lh_fct': None
                          }

    tree_distance_function = 'likelihood'
    if rec_parameter_dict['give_lh_fct'] is None:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['old_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['old_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['old_alpha']
    else:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['new_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['new_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['new_alpha']

    # data_path = os.path.join('data', 'omer_data')
    # data_path = os.path.join('data', '2211_alex_dataset', 'aligned_spacers_data',
    #                          'aligned_by_species_fc', 'dict_aligned_groups.pickle')
    data_path = os.path.join('data', '2211_alex_dataset', 'aligned_spacers_data',
                             # 'aligned_by_genus',
                             # 'aligned_by_genus_fixed_orientation',
                             # 'aligned_by_genus_pseudo_dna_fft',
                             # 'aligned_4_stdandfft_reindexed',
                             # 'aligned_5_stdandfft_reindexed_new_groups',
                             # 'aligned_7_stdandfft_reindexed_new_groups_combined_repeats_new_reindex',
                             'aligned_8_all_reindexed_data_groot',
                             'dict_aligned_groups.pickle')
    alignment_save_path = None
    selection_fct = None
    # (['crispr_type', 'minimum_array_number', 'selection not in groups'],
    #                      [['cas-i-e', 'cas-i-f', 'cas-ii-a', 'cas-i-c', 'cas-ii-c', 'cas-iii-a'], 2, ['group_561']],
    #                      )
    # selection_fct = [['selection of groups'], [['g_517']]]
    result_folder = os.path.join('0_result_folder', 'crispr_db')
    save_path = os.path.join(result_folder,
                             # 'determine_orientation_w_stdandfft_reindexed_new_groups_combined_repeats_new_reindex',
                             '0_all_data_on_groot',
                             )
    logfile_path = os.path.join(save_path, '0_logger.log')
    combine_non_unique_arrays = True
    tree_path = None
    tree_save_path = save_path
    hide_unobserved_spacers = False
    plot_order = False
    plot_tree = True
    significance_value = 0.05

    alignment_options_dict = {}

    experiments_run.run_pickled_data(rec_parameter_dict, data_path=data_path, save_path=save_path,
                                     plot_tree=plot_tree, logfile_path=logfile_path, do_show=False,
                                     combine_non_unique_arrays=combine_non_unique_arrays, tree_path=tree_path,
                                     tree_save_path=tree_save_path, hide_unobserved_spacers=hide_unobserved_spacers,
                                     selection_fct=selection_fct,
                                     plot_order=plot_order, significance_level=significance_value,
                                     finite_bl_at_root=finite_bl_at_root,
                                     alignment_options_dict=alignment_options_dict,
                                     group_by=group_by,
                                     extend_branches=extend_branches,
                                     tree_distance_function=tree_distance_function,
                                     tree_gain_rate=tree_gain_rate,
                                     tree_loss_rate=tree_loss_rate,
                                     tree_alpha=tree_alpha,
                                     determine_orientation=determine_orientation,
                                     orient_boundary=orient_boundary,
                                     )


def run_simulation_experiment():
    sim_as_rec = False
    logfile_path = os.path.join('testing_new_model', 'logger.log')
    save_path = 'testing_new_model'
    plot_tree = True
    parameter_dict = {'nb_trees': [10, 10],
                      'nb_leafs': [16, 16],
                      'gtr_gain_rate': [6, 6],
                      'deletion_model': ['independent', 'block'],
                      'lr_based_on_gr_avg_array_length': True,
                      'avg_array_length': [18, 18],
                      'deletion_model_parameter': [None, 1.5],
                      'event_sim': True,
                      'rnd_bl': True,
                      'tree_seed': 5,
                      'sim_seed': 5,
                      'model_name': 'g',
                      }
    experiments_run.run_simulation(parameter_dict, save_path, plot_tree=plot_tree,
                                   logfile_path=logfile_path,
                                   sim_as_rec=sim_as_rec)


def run_simulation_and_reconstruction_experiment():
    """
    Import from file sim_parameter_dict:
    deletion_model, import_path, repeat_runs, tree_sim_per_repeat_run (if trees_path is None), trees_path,
    tree_sim_seed, sim_seed, event_sim
    :return:
    """
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    # Old experiments
    # path = os.path.join('pictures', 'simulations', '15_blm_1_5_len_15_leafs_16_gr_6_new_lh_fct')
    # path = os.path.join('pictures', 'simulations', '16_blm_1_5_len_15_leafs_16_gr_6_old_lh_fct')
    # path = os.path.join('pictures', 'simulations', '17_blm_1_5_len_15_leafs_16_gr_6_new_lh_fct_sim_as_rec_w_unob')
    # path = os.path.join('pictures', 'simulations', '18_blm_1_5_len_15_leafs_16_gr_6_old_lh_fct_sim_as_rec_w_unob')
    # path = os.path.join('pictures', 'simulations', '19_blm_1_5_len_15_leafs_16_gr_6_new_lh_fct_sim_as_rec_wo_unob')
    # path = os.path.join('pictures', 'simulations', '20_blm_1_5_len_15_leafs_16_gr_6_old_lh_fct_sim_as_rec_wo_unob')
    # path = os.path.join('pictures', 'simulations', '21_blm_1_5_len_100_leafs_16_gr_6_new_lh_fct_sim_as_rec_w_unob')
    determine_orientation = True
    randomize_orientation = False
    orient_boundary = 10  # recommends reversion, or not,  or returns "uncertain".
    sim_simple_alignment = False
    finite_bl_at_root = True
    load_sim_from_pkl = False

    group_by = 'sim|loss model param.'
    do_show = False
    plot_order = False
    plot_tree = False
    result_folder = os.path.join('0_result_folder', 'sim_testing')
    # save_path = os.path.join(result_folder, '14_testing_likelihoodtree_omer_data_parameter')
    save_path = os.path.join(result_folder, 'test_data_from_pkl')
    logfile_path = os.path.join(save_path, '0_logger.log')
    rec_save_path = os.path.join(save_path, 'reconstruction')
    sim_save_path = os.path.join(result_folder, 'sim_data_test_data_from_pkl')
    # sim_save_path = os.path.join(result_folder, '0_simple_simulations')
    sim_logfile_path = os.path.join(save_path, '0_sim_logger.log')
    sim_as_rec = False
    hide_unobserved_spacers = False
    selection_fct = None
    significance_value = 0.05

    rec_parameter_dict = {'model': 'gtr', 'gain_rate': 1, 'loss_rate': 0.1,
                          # 'give_lh_fct': os.path.join('data', '0_lh', 'death_lh_up_to_100_lambdified.pkl'),
                          'give_lh_fct': None,
                          }
    use_sim_tree = True
    tree_distance_function = 'likelihood'
    if rec_parameter_dict['give_lh_fct'] is None:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['old_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['old_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['old_alpha']
    else:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['new_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['new_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['new_alpha']

    deletion_model_parameter = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.25, 2.5, 3.0, 4.0, 5.0]
    # deletion_model_parameter = [1.5]
    deletion_model = ['block'] * len(deletion_model_parameter)
    gtr_gain_rate = [6] * len(deletion_model_parameter)
    # gtr_gain_rate = [30] * len(deletion_model_parameter)
    nb_trees = [10] * len(deletion_model_parameter)
    # nb_trees = [10] * len(deletion_model_parameter)
    nb_leafs = [16] * len(deletion_model_parameter)
    # nb_leafs = [3] * len(deletion_model_parameter)
    avg_array_length = [18] * len(deletion_model_parameter)
    # avg_array_length = [8] * len(deletion_model_parameter)
    # avg_array_length = [100] * len(deletion_model_parameter)
    sim_parameter_dict = {'nb_trees': nb_trees,
                          'nb_leafs': nb_leafs,
                          'gtr_gain_rate': gtr_gain_rate,
                          'deletion_model': deletion_model,
                          'lr_based_on_gr_avg_array_length': True,
                          'avg_array_length': avg_array_length,
                          'deletion_model_parameter': deletion_model_parameter,
                          'event_sim': True,
                          'rnd_bl': True,
                          'tree_seed': 5,
                          'sim_seed': 5,
                          'model_name': 'g',
                          }

    import_folder = os.path.join('0_result_folder',
                                 'crispr_db',
                                 'determine_orientation_w_stdandfft_reindexed_new_groups_combined_repeats_new_reindex')
    # Note: blm_two_directional True + blm_rate_correction_last_spacer True leads to only corrected front spacer. I.e.
    # generally use either blm_two_directional or blm_rate_correction_last_spacer.
    sim_parameter_dict = {'deletion_model': 'block',
                          'event_sim': True,
                          'tree_seed': 5,
                          'sim_seed': 5,
                          'import_path': os.path.join(import_folder,
                                                      '0_final_oriented_protocol.pkl'),
                          'trees_path': os.path.join(import_folder, 'final_dict_nwk_trees.pickle'),
                          'repeat_runs': 2,
                          # 'blm_two_directional': True,
                          # 'blm_rate_correction_last_spacer': True,
                          }

    experiments_run.run_simulation_and_reconstruction(sim_parameter_dict, rec_parameter_dict, save_path=save_path,
                                                      rec_save_path=rec_save_path,
                                                      plot_tree=plot_tree, sim_save_path=sim_save_path,
                                                      sim_plot_tree=plot_tree, sim_logfile_path=sim_logfile_path,
                                                      logfile_path=logfile_path, do_show=do_show,
                                                      sim_as_rec=sim_as_rec,
                                                      hide_unobserved_spacers=hide_unobserved_spacers,
                                                      selection_fct=selection_fct, plot_order=plot_order,
                                                      significance_level=significance_value,
                                                      alignment_options=MAFFT_OPTIONS,
                                                      load_sim_from_pkl=load_sim_from_pkl,
                                                      finite_bl_at_root=finite_bl_at_root,
                                                      sim_simple_alignment=sim_simple_alignment,
                                                      group_by=group_by,
                                                      use_sim_tree=use_sim_tree,
                                                      tree_distance_function=tree_distance_function,
                                                      tree_gain_rate=tree_gain_rate,
                                                      tree_loss_rate=tree_loss_rate,
                                                      tree_alpha=tree_alpha,
                                                      determine_orientation=determine_orientation,
                                                      randomize_orientation=randomize_orientation,
                                                      orient_decision_boundary=orient_boundary
                                                      )


if __name__ == '__main__':
    pd.set_option('display.max_columns', None)
    pd.set_option('display.max_rows', None)
    sys.setrecursionlimit(20000)
    # args = config.load_config()
    # run_real_data()
    # run()
    # compare_msa()

    # run_simulation_experiment()
    # run_pickled_data()
    # run_crispr_db_data()
    # run_simulation_and_reconstruction_experiment()
