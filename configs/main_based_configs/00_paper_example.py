import experiments_run
import model.experiments_tools as experiments_tools
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


def run_crispr_db_data():
    """
    Backbone for simulations, reconstruction of whole dataset, excluding big groups. Naive likelihood function is used.
    :return:
    """
    determine_orientation = False
    finite_bl_at_root = True

    orient_boundary = 10
    extend_branches = False
    significance_value = 0.05
    group_by = None
    rec_parameter_dict = {'model': 'gtr', 'gain_rate': 10e-6, 'loss_rate': 10e-2,
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

    ############################################################################################################
    selection_fct = None
    # g_500 is example on poster!
    # selection_fct = [['maximum_unique_array_number'], [500]]
    # (['crispr_type', 'minimum_array_number', 'selection not in groups'],
    #                      [['cas-i-e', 'cas-i-f', 'cas-ii-a', 'cas-i-c', 'cas-ii-c', 'cas-iii-a'], 2, ['group_561']],
    #                      )
    ############################################################################################################
    result_folder = os.path.join('0_result_folder', '0_paper')
    data_path = os.path.join('data', '2_paper_example', '0_poster_group_3', '0_poster_group_3.fa')
    # tree_path = os.path.join('data', 'franz_test_tree_g13.pickle')
    tree_path = os.path.join('data', '2_paper_example', '0_poster_group_3', '0_toy_crispr_tree.pkl')
    core_genome_trees = False
    # tree_path = None
    save_path = os.path.join(result_folder,
                             # 'determine_orientation_w_stdandfft_reindexed_new_groups_combined_repeats_new_reindex',
                             '00_paper_example',
                             # '0_g_502_for_jan',
                             )
    tree_save_path = save_path
    logfile_path = os.path.join(save_path, '0_logger.log')
    ############################################################################################################
    alpha_bias_correction = True
    rho_bias_correction = True
    ############################################################################################################
    combine_non_unique_arrays = False
    hide_unobserved_spacers = False
    ############################################################################################################
    plot_order = 'graphviz'
    plot_tree = True
    ############################################################################################################
    alignment_options_dict = {}

    ############################################################################################################
    # Convert data to pickle
    sp = os.path.join('data', '2_paper_example', '0_poster_group_3', '0_poster_group_3.pickle')
    experiments_tools.load_align_single_fasta(data_path, os.path.join('data', '2_paper_example', '0_poster_group_3',
                                                                      'work_folder'),
                                              MAFFT_OPTIONS, group_name='s', save_path=sp)
    data_path = sp

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
                                     core_genome_trees=core_genome_trees,
                                     alpha_bias_correction=alpha_bias_correction,
                                     rho_bias_correction=rho_bias_correction,
                                     )


if __name__ == '__main__':
    run_crispr_db_data()
