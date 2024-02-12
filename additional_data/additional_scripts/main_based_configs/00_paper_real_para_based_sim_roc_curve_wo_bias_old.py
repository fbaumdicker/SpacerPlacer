import sys

import experiments_run
import os

############################################################################################################
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

############################################################################################################
# Update those? What should be the basis?
TREE_GENERATION_PARAMETERS = {'old_gain_rate': 0.9724,  # based on omer data (without boring data; old and new lh fct)!
                              'old_loss_rate': 0.1765,
                              'old_alpha': 5.5092,
                              'new_gain_rate': 0.6108,
                              'new_loss_rate': 0.1830,
                              'new_alpha': 3.3377,
                              }


def determine_orientation_threshold():
    """
    Determine sensible orientation threshold based on basic dataset.
    :return:
    """
    finite_bl_at_root = True
    determine_orientation = True
    randomize_orientation = False
    orient_boundary = 0
    plot_tree = False
    plot_order = False
    do_show = False
    tree_distance_function = 'likelihood'
    # lrt test significance value
    significance_value = 0.05

    extend_branches = 10e-8

    # selection_fct = None
    selection_fct = [['empty_arrays'], [None]]
    ############################################################################################################
    # paths
    # result_folder = os.path.join('0_result_folder', '0_paper', '0_rec_paper', '0_determine_orientation_threshold')
    result_folder = os.path.join('0_result_folder', '0_paper')
    experiment_name = '00_paper_real_para_based_sim_roc_curve_wo_bias_old'
    save_path = os.path.join(result_folder, experiment_name)
    logfile_path = os.path.join(save_path, '0_logfile.log')
    sim_save_path = os.path.join(save_path, '1_simulations')
    rec_save_path = os.path.join(save_path, '2_reconstruction')
    tree_save_path = os.path.join(save_path, '3_trees')
    sim_logfile_path = os.path.join(sim_save_path, '0_sim_logfile.log')

    import_folder = os.path.join('0_result_folder', '0_paper', '1_excluding_big_dataset_n_lh_fct_cg_trees')
    # import_folder = os.path.join('0_result_folder', '0_paper',
    #                              '00_paper_real_data_optimal_estimation_alt')
    sim_parameter_dict = {'deletion_model': 'block',
                          'event_sim': True,
                          'tree_seed': 500,  # previously 5
                          'sim_seed': 500,  # previously 5
                          'import_path': os.path.join(import_folder, '0_final_oriented_protocol.pkl'),
                          'trees_path': os.path.join(import_folder, 'final_dict_nwk_trees.pickle'),
                          'repeat_runs': 10,
                          'blm_rate_correction_last_spacer': False,
                          'blm_two_directional': True,
                          }
    ############################################################################################################
    # bias correction
    alpha_bias_correction = True
    rho_bias_correction = True
    ############################################################################################################
    # Simulation already available?
    # load_sim_from_pkl = True
    load_sim_from_pkl = False
    ############################################################################################################
    # Use simulated trees or estimate...
    # use_sim_tree = False
    use_sim_tree = True
    ############################################################################################################
    # Simulation as reconstruction ....
    # sim_as_rec = True
    sim_as_rec = False
    ############################################################################################################
    # hide unobserved spacers:
    hide_unobserved_spacers = True
    # hide_unobserved_spacers = False

    ############################################################################################################
    # simulation simple alignment:
    # sim_simple_alignment = False
    sim_simple_alignment = True

    ############################################################################################################
    # BDM Naive vs. BDM ode based:
    ode_based_lh_path = os.path.join('data', '0_lh', '230329_death_lh_up_to_100_groot.pickle')
    rec_parameter_dict = {'model': 'gtr', 'gain_rate': 1e-5, 'loss_rate': 1e-1,
                          'give_lh_fct': None,
                          # 'give_lh_fct': ode_based_lh_path,
                          }
    if rec_parameter_dict['give_lh_fct'] is None:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['old_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['old_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['old_alpha']
    else:
        tree_gain_rate = TREE_GENERATION_PARAMETERS['new_gain_rate']
        tree_loss_rate = TREE_GENERATION_PARAMETERS['new_loss_rate']
        tree_alpha = TREE_GENERATION_PARAMETERS['new_alpha']

    ############################################################################################################
    # group_by for estimates over large parts of dataset
    group_by = 'sim|repeat run'

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
                                                      orient_decision_boundary=orient_boundary,
                                                      alpha_bias_correction=alpha_bias_correction,
                                                      rho_bias_correction=rho_bias_correction,
                                                      extend_branches=extend_branches,
                                                      )


if __name__ == '__main__':
    sys.setrecursionlimit(10000)
    determine_orientation_threshold()
