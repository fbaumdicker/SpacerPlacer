import os


def compose_summary_dict(df_results, supplied_parameters):
    estimated_parameters_keys = {'deletion_rate_bdm', 'alpha_bdm', 'per_spacer_deletion_rate_bdm', 'insertion_rate_bdm',
                                 'deletion_rate_idm',
                                 'insertion_rate_idm',
                                 }
    idm_vs_bdm_test_keys = {'ln_lh_idm', 'ln_lh_bdm', 'test statistic (-2*(ln_lh_ratio)', 'test result',
                            'Deletion model preferred by LRT',
                            'chi2_quantile', 'significance_level'}
    reconstruction_results_keys = {'nb of spacers in alignment', 'nb of unique spacers',
                                   'nb of reconstructed insertions',
                                   'nb of reconstructed deletions',
                                   'nb of reconstructed dup/rearr candidates',
                                   'nb of reconstructed rearrangements',
                                   'nb of reconstructed duplications',
                                   'nb of reconstructed reacquisitions',
                                   'nb of reconstructed independent gains',
                                   'nb of reconstructed other dup. events',
                                   'nb of reconstructed ind. acquisitions, not duplicates'}
    orientation_results_keys = {'reversed_ln_lh_idm', 'reversed_ln_lh_bdm', 'ln_lh_bdm - reversed_ln_lh_bdm',
                                'predicted orientation'}
    tree_information_keys = {'tree height', 'tree length', 'min/max tree branch lengths',
                             'nb of leafs (before combining non-uniques)',
                             'nb of leafs (after combining non-uniques)',
                             }
    general_information_keys = {'nb of unique spacer arrays',
                                'run time'}

    dict_results = df_results.to_dict(orient='index')
    ordered_results = dict()

    for key, value in dict_results.items():
        est_par = {k: v for k, v in value.items() if k in estimated_parameters_keys}
        idm_vs_bdm_test = {k: v for k, v in value.items() if k in idm_vs_bdm_test_keys}
        rec_res = {k: v for k, v in value.items() if k in reconstruction_results_keys}
        orient_res = {k: v for k, v in value.items() if k in orientation_results_keys}
        tree_info = {k: v for k, v in value.items() if k in tree_information_keys}
        gen_info = {k: v for k, v in value.items() if k in general_information_keys}
        ordered_results[key] = {'estimated_parameters': est_par, 'idm_vs_bdm_test': idm_vs_bdm_test,
                                'reconstruction_results': rec_res, 'orientation_results': orient_res,
                                'tree_information': tree_info, 'general_information': gen_info}

    dict_summary = {'input_path': supplied_parameters.pop('input_path'),
                    'output_path': supplied_parameters.pop('output_path'), 'supplied_parameters': supplied_parameters,
                    'ind_groups': ordered_results}

    return dict_summary


def write_summary(dict_summary, out_path):
    determine_orientation = dict_summary['supplied_parameters']['determine_orientation']
    combine_non_uniques = dict_summary['supplied_parameters']['combine_non_unique_arrays']

    results_path = os.path.join(dict_summary['output_path'], '0_results.csv')
    detailed_results_path = os.path.join(dict_summary['output_path'], 'detailed_results')
    additional_data_path = os.path.join(dict_summary['output_path'], 'additional_data')
    forward_plots_path = os.path.join(dict_summary['output_path'], '0_forward')
    reversed_plots_path = os.path.join(dict_summary['output_path'], '0_reversed')

    with open(out_path, 'w') as file:
        file.write(f'Reconstruction was successful!\n')
        file.write(f'Input and Output locations:\n')
        file.write(f'\tInput path: {dict_summary["input_path"]}\n')
        file.write(f'\tOutput path: {dict_summary["output_path"]}\n')
        file.write('\n')

        file.write('More details can be found in the following folders/files:\n')
        file.write(f'\tResults in tabular format: {results_path}\n')
        file.write(f'\tMore detailed results: {detailed_results_path}\n')
        file.write(f'\tAdditional data such as spacer fasta files and trees: {additional_data_path}\n')
        if determine_orientation:
            file.write(f'\tPlots and detailed results for forward orientation: {forward_plots_path}\n')
            file.write(f'\tPlots and detailed results for reversed orientation: {reversed_plots_path}\n')
        else:
            file.write(f'\tPlots: {forward_plots_path}\n')
        file.write('\n')
        file.write('Used parameters:\n')
        for key, value in dict_summary['supplied_parameters'].items():
            file.write(f'\t{key}: {value}\n')
        file.write('\n')

        file.write('Result summary for each reconstructed group:\n')
        file.write('\n')
        for group, dict_sum in dict_summary['ind_groups'].items():
            file.write(f'{group}:\n')
            file.write(f'\n')
            file.write(f'\tEstimated parameters:\n')
            for key, value in dict_sum['estimated_parameters'].items():
                file.write(f'\t{key}: {value}\n')
            file.write(f'\n')
            file.write(f'\tLikelihood ratio test:\n')
            for key, value in dict_sum['idm_vs_bdm_test'].items():
                file.write(f'\t{key}: {value}\n')
            file.write(f'\n')
            file.write(f'\tReconstruction details:\n')
            for key, value in dict_sum['reconstruction_results'].items():
                file.write(f'\t{key}: {value}\n')
            if determine_orientation:
                file.write(f'\n')
                file.write(f'\tOrientation prediction:\n')
                for key, value in dict_sum['orientation_results'].items():
                    # if key == 'recommend reversing array':
                    #     value = 'Reverse' if value is True else 'Forward'
                    #     file.write(f'Predicted orientation: {value}\n')
                    file.write(f'\t{key}: {value}\n')
            file.write(f'\n')

            if not combine_non_uniques:
                _ = dict_sum['tree_information'].pop('nb of leafs (after combining non-uniques)')
                nb_leafs = dict_sum['tree_information'].pop('nb of leafs (before combining non-uniques)')
                dict_sum['tree_information']['nb of leafs'] = nb_leafs

            file.write(f'\tTree information:\n')
            for key, value in dict_sum['tree_information'].items():
                file.write(f'\t{key}: {value}\n')
            file.write(f'\n')
            file.write(f'\tGeneral information:\n')
            for key, value in dict_sum['general_information'].items():
                file.write(f'\t{key}: {value}\n')
        return
