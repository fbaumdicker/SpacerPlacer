import argparse
import os
import json

import experiments_run
from input_parser import InputParser
from model.helpers.misc import create_logger
from model.summary import compose_summary_dict, write_summary
from experiments_run import run_pickled_data, run_multiple_groups
import model.helpers.import_data as import_data


class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


parser = ArgumentParser(prog='spacerplacer',
                        description='SpacerPlacer: An Ancestral Reconstruction Algorithm for CRISPR Arrays',
                        fromfile_prefix_chars='@',
                        )
############################################################################################################
# General:
# parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--verbosity', type=int, default=2, choices=[0, 1, 2],
                    help='Verbosity level. 0: no output, 1: minimal output, 2: maximal output.')
parser.add_argument('--seed', type=int, default=2357, help='Seed for the random number generator '
                                                           'only impact on MAFFT.')
# parser.add_argument('--group_by_metadata', type=str, default=None,
#                     help='Not supported.')
############################################################################################################
# Input/Output locations + format:
parser.add_argument('input_path', type=str,
                    help='Path to input file or folder containing input files. The required input format is described '
                         'in the readme.')
parser.add_argument('output_path', type=str,
                    help='Path to output folder. The output folder will be created if it does not exist. '
                         'The output folder will contain files as described in the readme.')
parser.add_argument('--tree_path', type=str, default=None,
                    help='Path to tree json file or folder with newick files. '
                         'If none is provided, trees are estimated by SpacerPlacer. '
                         'The trees can be given in a newick format or as a dictionary in a json file '
                         '(such a json file is returned by SpacerPlacer).')
parser.add_argument('-it', '--input_type', type=str, choices=['spacer_fasta', 'pickled', 'ccf',
                                                              'crisprcasfinder'],
                    default='spacer_fasta',  # think about this name
                    help='Determines the input type, i.e. either already preprocessed fasta style spacer arrays, '
                         'a pickled file with CRISPR group(s) (our own data structure) '
                         'or data extracted from CRISPRCasFinder or CRISPRCasdb.')
# Have to decide on this at some point:
# parser.add_argument('--metadata_arrays', type=str, default=None,
#                     help='Path to metadata file. The metadata file should be a csv or dictionary')
# parser.add_argument('--metadata_spacers', type=str, default=None,
#                     help='Path to metadata file. The metadata file should be a csv or dictionary')
# parser.add_argument('-mg', '--multiple_groups', action='store_true', default=False,
#                     help='If given, flags that the analysis of multiple groups is desired.')

############################################################################################################
# Reconstruction:
parser.add_argument('--rec_model', type=str, choices=['gtr', ],
                    default='gtr', help='Determines the reconstruction model. Currently redundant, '
                                        'as only "gtr" is implemented.')
parser.add_argument('--insertion_rate', type=float, default=1e-5,
                    help='Insertion rate of the reconstruction model. Generally should be chosen to be small '
                         'compared to the deletion rate. Otherwise may lead to high number of independent '
                         'acquisitions and worse reconstructions. The ratio between insertion and deletion rate '
                         'can be tuned in case of reconstructions with excessive independent acquisitions or excessive '
                         'accumulations of insertions at the root. '
                         'Default insertion and deletion rates are provided and '
                         'work well for the tested datasets.')
parser.add_argument('--deletion_rate', type=float, default=1e-1,
                    help='Deletion rate of the reconstruction model. Generally should be chosen to be large '
                         'compared to the insertion rate. Otherwise may lead to high number of independent '
                         'acquisitions and worse reconstructions. The ratio between insertion and deletion rate '
                         'can be tuned in case of reconstructions with excessive independent acquisitions or excessive '
                         'accumulations of insertions at the root. '
                         'Default insertion and deletion rates are provided and '
                         'work well for the tested datasets.')
# Don't know if this is needed:
parser.add_argument('--extend_branches', type=float, default=0.00001,
                    help='Extends branches of the tree by the given value. This is useful, if the tree is not '
                         'well resolved to allow placement of events on '
                         'small (or length=0) branches.')
# parser.add_argument('--stationary_dist_at_root', action='store_true', default=False,
#                     help="If given, the likelihood at the root is computed based on the stationary distribution. "
#                          "Otherwise the maximum of the child's branch lengths and the standard likelihood "
#                          "computation is used. The stationary distribution might lead to undesired results, "
#                          "that are not parsimonious. We do not recommend to use this flag.")
# ########################################################################################################### Tree
# tree estimation:
parser.add_argument('--tree_distance_function', type=str, choices=['likelihood'],
                    default='likelihood', help='Determines the distance function used for the tree estimation '
                                               'with UPGMA. Currently redundant, as only "likelihood" is implemented.')
parser.add_argument('--tree_construction_method', type=str, choices=['upgma', 'nj'],
                    default='upgma', help='Determines the tree construction method used for the tree estimation. '
                                          'Currently UPGMA (upgma) and neighbor joining (nj) are implemented.')
parser.add_argument('--tree_insertion_rate', type=float, default=None,
                    help='The user can provide their own insertion rate for the tree estimation based on the '
                         'block deletion likelihood function. '
                         'In the default case, well performing parameters are chosen '
                         'depending on the used likelihood function ("simple", "ode_based").')
parser.add_argument('--tree_deletion_rate', type=float, default=None,
                    help='The user can provide their own deletion rate for the tree estimation based on the '
                         'block deletion likelihood function. '
                         'In the default case, well performing parameters are chosen '
                         'depending on the used likelihood function ("simple", "ode_based").'
                    )
parser.add_argument('--tree_alpha', type=float, default=None,
                    help='The user can provide their own deletion rate for the tree estimation based on the '
                         'block deletion likelihood function. '
                         'In the default case, well performing parameters are chosen '
                         'depending on the used likelihood function ("simple", "ode_based").'
                    )
parser.add_argument('--tree_lh_fct', type=str, choices=['simple', 'ode_based'], default='simple',
                    help='Determines the likelihood function used for the tree estimation with UPGMA. '
                         'Details about the likelihood functions can be found in the paper. The default "simple" has '
                         'been found to perform better in simulations.')
parser.add_argument('--combine_non_unique_arrays', action='store_true',
                    help='If given, arrays with the exactly the same spacer content are combined before the '
                         'tree estimation and reconstruction. This might be helpful reduce clutter '
                         '(especially in visualization).'
                         'Arrays with the same spacer content can not be separated in the tree '
                         'estimation anyway. '
                         'The default is NOT to combine arrays with the same spacer content. '
                         'If arrays are combined, the combined array names are saved in detailed results csv under '
                         'column "combined array names".')

############################################################################################################
# Visualization:
parser.add_argument('--no_plot_reconstruction', action='store_true',
                    help='If given, the reconstruction is not plotted.')
parser.add_argument('--no_plot_order_graph', action='store_true',
                    help='If given, the Partial Spacer Insertion Order is not plotted.')
# Does this work?
parser.add_argument('--do_show', action='store_true',
                    help='If given, the plots are shown directly.')

############################################################################################################
# Evaluation (LRT parameters):
parser.add_argument('--lh_fct', type=str, choices=['simple', 'ode_based'], default='simple',
                    help='Determines the likelihood function used for the Block deletion model estimates and '
                         'the likelihood ratio test. "simple" allows bias '
                         'corrections in rho and alpha and was found to perform better in simulations. '
                         '"ode_based" is the less performant, but, in theory, more precise likelihood function.')
parser.add_argument('--no_alpha_bias_correction', action='store_true',
                    help='If given, omits the alpha bias correction. Best performance is achieved with '
                         'the bias correction. Bias correction is only possible, '
                         'if simple likelihood function is used.')
parser.add_argument('--no_rho_bias_correction', action='store_true',
                    help='If given, omits the rho bias correction. Best performance is achieved with '
                         'the bias correction. Bias correction is only possible, '
                         'if simple likelihood function is used.')
parser.add_argument('--significance_level', type=float, default=0.05,
                    help='Determines the significance level for the likelihood ratio test between the Independent '
                         'Deletion Model and the Block Deletion Model. '
                         'The test statistic (-2*ln(lh_idm/ln_bdm) = 2*ln(lh_bdm - ln_idm) is chi-squared '
                         "distributed with one single degree of freedom by Wilks' Theorem.")

############################################################################################################
# Orientation determination:
# Do we want to always determine the orientation (for CRISPRCasdb data)? Or make it dependent on the input structure?
parser.add_argument('--determine_orientation', action='store_true',
                    help='If given, both a forward and a reverse '
                         'reconstruction is made and the orientation of the arrays is determined by SpacerPlacer. '
                         'Otherwise, the orientation is '
                         'accepted as provided. Note, if the tree is estimated by SpacerPlacer, '
                         'a separate tree ist estimated for the reversed array, i.e. the trees between forward and '
                         'reverse orientation might differ.'
                    )
parser.add_argument('--orientation_decision_boundary', type=float, default=5,
                    help='Determines the decision boundary for the orientation decision. If the modulus of the '
                         'difference between the likelihoods of the two orientations is smaller than '
                         'the decision boundary, the orientation is set to "not determined" (ND). '
                         'Otherwise, the orientation '
                         'is set to the orientation with the higher likelihood. The decision boundary is given in '
                         'logscale.'
                         'The default value %(default)s was determined empirically and works well for the tested '
                         'datasets.')
############################################################################################################
# Additional data:
parser.add_argument('--save_reconstructed_events', action='store_true',
                    help='If given, the tree and reconstructed events along the tree are saved in '
                         '"reconstructed_events". On the basis of this data, the reconstruction is '
                         'visualized and the events can be analyzed in more detail. Currently only works, if the '
                         'the reconstruction is visualized.')
############################################################################################################
# Preprocessing:
# parser.add_argument('--clustering', action='store_true',
#                     help='If given, the input data is clustered in subgroups before the reconstruction.'
#                          ' The clustering process is described in the paper.')
# parser.add_argument('--group_by', type=str, choices=['crispr type',
#                                                      'minimum array number', 'selection not in groups',])

args = parser.parse_args()

if args.verbosity == 0:
    logging_level = 'WARNING'
elif args.verbosity == 1:
    logging_level = 'INFO'
else:
    logging_level = 'DEBUG'
logfile_path = os.path.join(args.output_path, '0_logger.log')

if not os.path.exists(args.output_path):
    os.makedirs(args.output_path)
logger = create_logger('spacerplacer', args.verbosity, outfile=logfile_path)

logger.info(f'Starting SpacerPlacer on Input path: {args.input_path}')
# start with function for single group -> multiple groups
logger.debug(f'Input type: {args.input_type}')

rec_parameter_dict = {'model': args.rec_model, 'gain_rate': args.insertion_rate, 'loss_rate': args.deletion_rate,
                      }

logger.debug(vars(args))

tree_path = None
extend_branches = args.extend_branches

if args.input_type == 'pickled':
    dict_crispr_groups = run_pickled_data(rec_parameter_dict,
                                          lh_fct=args.lh_fct,
                                          data_path=args.input_path,
                                          save_path=args.output_path,
                                          plot_tree=not args.no_plot_reconstruction,
                                          logfile_path=os.path.join(args.output_path, '0_logger.log'),
                                          plot_order=not args.no_plot_order_graph,
                                          combine_non_unique_arrays=False,
                                          tree_path=args.tree_path,
                                          tree_save_path=args.output_path,
                                          hide_unobserved_spacers=False,
                                          selection_fct=None,
                                          significance_level=args.significance_level,
                                          finite_bl_at_root=True,
                                          alignment_options_dict=None,
                                          group_by=None,
                                          extend_branches=extend_branches,
                                          tree_distance_function=args.tree_lh_fct,
                                          tree_gain_rate=args.tree_insertion_rate,
                                          tree_loss_rate=args.tree_deletion_rate,
                                          tree_alpha=args.tree_alpha,
                                          determine_orientation=args.determine_orientation,
                                          orient_boundary=args.orientation_decision_boundary,
                                          alpha_bias_correction=not args.no_alpha_bias_correction,
                                          rho_bias_correction=not args.no_rho_bias_correction,
                                          core_genome_trees=True if args.tree_path is not None else False,
                                          save_reconstructed_events=args.save_reconstructed_events, )
elif args.input_type in ['ccf', 'crisprcasfinder', 'spacer_fasta']:
    if args.input_type == 'spacer_fasta':
        if os.path.splitext(args.input_path)[-1] in {'.fa', '.fasta', '.fna'}:
            ls_path_to_spacer_fasta = [args.input_path]
        else:
            ls_path_to_spacer_fasta = [os.path.join(args.input_path, group) for group in os.listdir(args.input_path)
                                       if group.split('.')[-1] in {'fa', 'fasta', 'fna'}]

        if not ls_path_to_spacer_fasta:
            logger.error(f'No fasta file found in {args.input_path}.')
            raise ValueError(f'No fasta file found in {args.input_path}.')

    else:
        ls_path_to_spacer_fasta = []
        path_to_spacer_fasta_folder = os.path.join(args.output_path, 'additional_data', 'spacer_fasta')
        if not os.path.exists(path_to_spacer_fasta_folder):
            os.makedirs(path_to_spacer_fasta_folder)
        for group in os.listdir(args.input_path):
            path_to_spacer_fasta = os.path.join(path_to_spacer_fasta_folder, group + '.fa')
            path_to_n_seq_file = os.path.join(path_to_spacer_fasta_folder,
                                              group + '_spacer_name_to_seq.fa')
            ccf_parser = InputParser(os.path.join(args.input_path, group), path_to_spacer_fasta,
                                     spacer_number_to_seq_file=path_to_n_seq_file)
            ls_path_to_spacer_fasta.append(path_to_spacer_fasta)

    path_to_tree = os.path.join(args.output_path, 'additional_data')
    tree_path = os.path.join(args.output_path, 'additional_data', 'dict_nwk_trees.json')
    if args.tree_path is None:
        tree_path = None
    elif os.path.isfile(args.tree_path):
        logger.info(f'Using tree file {args.tree_path}.')
        if os.path.splitext(args.tree_path)[1] in ['.json', '.pkl', '.pickle']:
            tree_path = args.tree_path
        elif len(ls_path_to_spacer_fasta) == 1:
            dict_trees = {os.path.splitext(os.path.basename(ls_path_to_spacer_fasta[0]))[0]:
                              import_data.load_single_tree(args.tree_path).format(fmt='newick')}
            if not os.path.exists(os.path.dirname(path_to_tree)):
                os.makedirs(os.path.dirname(path_to_tree))
            json.dump(dict_trees, open(tree_path, 'w'))
    else:
        logger.info(f'Using tree folder {args.tree_path}.')
        dict_trees = {}
        list_dir_os_path = os.listdir(args.tree_path)
        list_dir_os_path_no_ext = set([os.path.splitext(o)[0] for o in os.listdir(args.tree_path)])
        for group in os.listdir(args.input_path):
            group_no_ext = os.path.splitext(group)[0]
            tp = None
            for file in list_dir_os_path:
                if file.startswith(group_no_ext) and not file.endswith(('.fa', '.fasta', '.fna')):
                    tp = os.path.join(args.tree_path, file)
                    break

            if tp is None:
                logger.error(f'No tree found for group {group}. '
                             f'File would be found in folder {args.tree_path} '
                             f'with filename {group_no_ext} (+ extension).')
                raise ValueError(f'No tree found for group {group}. '
                                 f'File would be found in folder {args.tree_path} '
                                 f'with filename {group_no_ext} (+ extension).')


            dict_trees[group_no_ext] = import_data.load_single_tree(tp).format(fmt='newick')
        if not os.path.exists(os.path.dirname(path_to_tree)):
            os.makedirs(os.path.dirname(path_to_tree))
        json.dump(dict_trees, open(tree_path, 'w'))

    df_results_wo_details = run_multiple_groups(ls_path_to_spacer_fasta, args.output_path,
                                                rec_parameter_dict,
                                                logger=logger,
                                                lh_fct=args.lh_fct,
                                                plot_tree=not args.no_plot_reconstruction,
                                                do_show=args.do_show,
                                                determine_orientation=args.determine_orientation,
                                                orientation_decision_boundary=args.orientation_decision_boundary,
                                                tree_path=tree_path,
                                                plot_order=not args.no_plot_order_graph,
                                                significance_level=args.significance_level,
                                                extend_branches=extend_branches,
                                                tree_distance_fct=args.tree_distance_function,
                                                tree_construction_method=args.tree_construction_method,
                                                tree_lh_fct=args.tree_lh_fct,
                                                tree_insertion_rate=args.tree_insertion_rate,
                                                tree_deletion_rate=args.tree_deletion_rate,
                                                tree_alpha=args.tree_alpha,
                                                alpha_bias_correction=not args.no_alpha_bias_correction,
                                                rho_bias_correction=not args.no_rho_bias_correction,
                                                combine_non_unique_arrays=args.combine_non_unique_arrays,
                                                seed=args.seed,
                                                save_reconstructed_events=args.save_reconstructed_events,
                                                )
    summary_dict = compose_summary_dict(df_results_wo_details, dict(vars(args)))
    write_summary(summary_dict, os.path.join(args.output_path, 'summary.txt'))

else:
    raise logger.error(f'Input type {args.input_type} not recognized.')
