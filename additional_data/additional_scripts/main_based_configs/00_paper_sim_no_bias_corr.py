import numpy as np

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


############################################################################################################
# LRT quality test:
def run_lrt_quality_test():
    finite_bl_at_root = True
    determine_orientation = False
    randomize_orientation = False
    orient_boundary = 0
    plot_tree = False
    plot_order = False
    do_show = False
    tree_distance_function = 'likelihood'
    # lrt test significance value
    significance_value = 0.05

    selection_fct = None
    # selection_fct = [['empty_arrays'], [None]]
    ############################################################################################################
    # paths
    result_folder = os.path.join('0_result_folder', '0_paper')
    experiment_name = '00_paper_sim_no_bias_corr'
    save_path = os.path.join(result_folder, experiment_name)
    logfile_path = os.path.join(save_path, '0_logfile.log')
    sim_save_path = os.path.join(save_path, '1_simulations')
    sim_logfile_path = os.path.join(sim_save_path, '0_sim_logfile.log')
    rec_save_path = os.path.join(save_path, '2_reconstruction')
    tree_save_path = os.path.join(save_path, '3_trees')
    ############################################################################################################
    # bias correction variables
    alpha_bias_correction = False
    rho_bias_correction = False
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
    # group_by for estimates over large parts of dataset
    group_by = 'sim|loss model param.'

    ls_alpha = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8,
                2.9, 3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8,
                3.9, 4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7,
                4.8, 4.9, 5]
    dataset_avg_array_length = 22  # 18  # last dataset had ~22.38, median ~15.625
    dataset_gain_rate = 6
    sim_parameter_dict = {
        'nb_trees': [1000] * len(ls_alpha),
        'deletion_model': ['block'] * len(ls_alpha),
        'nb_leafs': [16] * len(ls_alpha),  # ?????
        'deletion_model_parameter': ls_alpha,
        'gtr_gain_rate': [dataset_gain_rate] * len(ls_alpha),
        'avg_array_length': [dataset_avg_array_length] * len(ls_alpha),
        'event_sim': True,
        'rnd_bl': True,
        'tree_seed': 5,
        'sim_seed': 5,
        'model_name': 'g',
        'lr_based_on_gr_avg_array_length': True,
        'blm_rate_correction_last_spacer': False,
        'blm_two_directional': True,
    }

    ############################################################################################################
    # BDM Naive vs. BDM ode based:
    ode_based_lh_path = os.path.join('data', '0_lh', '230329_death_lh_up_to_100_groot.pickle')
    rec_parameter_dict = {'model': 'gtr', 'gain_rate': 1, 'loss_rate': 0.1,
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
                                                      )


if __name__ == '__main__':
    run_lrt_quality_test()
