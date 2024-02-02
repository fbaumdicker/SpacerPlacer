import subprocess
import os
import pickle

from model.data_classes.advanced_tree import AdvancedTree
from model.data_classes.crisprdata import CRISPRGroup
from model.helpers import import_data

MAFFT_OPTIONS = [  # '--genafpair',
    '--localpair',
    # '--sp',
    # '--globalpair',
    # '--hybridpair',
    # '--leavegappyregion',
    '--maxiterate', '1000',
    '--randomseed', '2357',
    '--quiet',
    '--thread', '-1',  # multithreading -1 -> automatically choose
    # '--threadit', '0',  # threaded iterative alignment has some randomness (this stops this)
]


def root_tree(tree_path, root_name):
    if os.path.exists(tree_path):
        tree = import_data.load_single_tree(tree_path)
    else:
        return None
    tree.root_at_midpoint()
    tree.root.name = root_name
    return tree


def create_cas_sequence_files(dict_chosen_cas, save_path, group_name='g'):
    """
    Cas sequences are lists, their labeling is given by cas_info.
    :param dict_cas: Note, sequences are Bio.Seq
    :param save_path:
    :return:
    """
    dict_s_p = dict()
    for cas_name, dict_cas in dict_chosen_cas.items():
        s_p = os.path.join(save_path, group_name + '_' + cas_name)
        dict_s_p[cas_name] = s_p
        with open(s_p, 'w') as f:
            for ind_n, seq in dict_cas.items():
                f.write('>' + ind_n)
                f.write('\n')
                f.write(str(seq))
                f.write('\n')
    return dict_s_p


def create_cas_sequence_file_2(crispr_group, save_path):
    """
    Cas sequences are lists, their labeling is given by cas_info.
    :param crispr_group:
    :param save_path:
    :return:
    """
    dict_cas_info = crispr_group.all_cas_types()
    dict_cas_type = crispr_group.get_cas_type()
    for name, info in dict_cas_info.items():
        print(name, dict_cas_type[name], info)
    print()
    dict_cas_seq = crispr_group.pop_cas_sequences()
    with open(save_path, 'w') as f:
        for name, cas_seq in dict_cas_seq.items():
            f.write('>' + name)
            f.write('\n')
            f.write(''.join(cas_seq))
            f.write('\n')
    return


def combine_mafft_sequences(dict_input_paths, save_path):
    dict_aligned_cas_genes_by_name = dict()
    for cas_name, path in dict_input_paths.items():
        with open(path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if line[0] == '>':
                    name = line[1:].strip('\n')
                else:
                    seq = line.strip('\n')
                    # seq = seq.replace('--', '-')
                    if name in dict_aligned_cas_genes_by_name:
                        dict_aligned_cas_genes_by_name[name] += '-' + seq
                    else:
                        dict_aligned_cas_genes_by_name[name] = seq
    with open(save_path, 'w') as f:
        for ind_n, cas_genes in dict_aligned_cas_genes_by_name.items():
            f.write('>' + ind_n)
            f.write('\n')
            f.write(cas_genes)
            f.write('\n')


def run_mafft_for_cas_sequences(dict_input_paths, save_path_aligned_by_name, mafft_options=None):
    """
    L-INS-i is supposed to be the best performing mafft model (in general). Speed is probably not important.
    :param input_path:
    :param save_path:
    :return:
    """
    mafft_options = MAFFT_OPTIONS if mafft_options is None else mafft_options
    dict_s_p = dict()
    for cas_name, input_path in dict_input_paths.items():
        s_p = input_path + '_mafft_aligned'
        dict_s_p[cas_name] = s_p
        with open(s_p, 'w') as outfile:
            print('Calling ', ' '.join(['mafft'] + [input_path]))
            subprocess.call(['mafft'] + mafft_options + [input_path], stdout=outfile)
    combine_mafft_sequences(dict_s_p, save_path_aligned_by_name)
    return dict_s_p


def run_iqtree(input_path, options=None):
    """
    Collins Paper used ModelFinder to suggest the model, so why not.... Do I want additional bootstrapping ('-bb')?
    :param options:
    :param input_path:
    :param save_path:
    :return:
    """
    if options is None:
        options = ['-s', input_path, '-bb', '1000', '-nt', 'AUTO']
    else:
        options = ['-s', input_path] + options
    print('Calling ', ' '.join(options))
    subprocess.call(['iqtree'] + options)


def create_cas_tree(crispr_group, work_path, save_path):
    create_cas_sequence_files(crispr_group, os.path.join(work_path, crispr_group.repeat + '.txt'))

    # run_mafft_for_cas_sequences(os.path.join(work_path, crispr_group.repeat + '.txt'),
    #                                      os.path.join(work_path, crispr_group.repeat + '_aligned.txt'))

    # run_iqtree(os.path.join(work_path, crispr_group.repeat + '_aligned.txt'),
    #            os.path.join(save_path, crispr_group.repeat + '.nwk'))
    return


if __name__ == '__main__':
    options = ['-nt', 'AUTO'] + ['-m', 'GTR+F+R9']  #  'AUTO' '-bb', '1000',
    coralignments_path = os.path.join('data', 'aligned_8_corealignments')
    dict_trees = {}
    all_fa_files = [os.path.join(coralignments_path, f) for f in os.listdir(coralignments_path) if f.endswith('.fa')]

    for fa_file in all_fa_files:
        with open(fa_file, 'r') as f:
            lines = f.readlines()
            count = 0
            for line in lines:
                if line[0] == '>':
                    count += 1
        if count < 4:
            continue
        if not os.path.exists(fa_file + '.treefile'):
            run_iqtree(fa_file, options=options)

        tree = import_data.load_single_tree(fa_file + '.treefile')
        tree = AdvancedTree(tree, False)
        # print(tree)

        group_name = fa_file.split('/')[-1].split('.')[0]
        group_name = '_'.join(group_name.split('_')[1:])
        dict_trees[group_name] = tree.format(format='newick')

        with open(os.path.join('data', 'aligned_8_corealignment_trees.pickle'), 'wb') as f:
            pickle.dump(dict_trees, f)
    # data_folder = os.path.join('data', '2211_alex_dataset')
    # save_path = os.path.join(data_folder, 'iqtree_reconstructed_trees')
    # work_path = os.path.join(data_folder, 'iqtree_work')
    # if not os.path.exists(save_path):
    #     os.makedirs(save_path)
    # if not os.path.exists(work_path):
    #     os.makedirs(work_path)
    # raw_data_path = os.path.join(data_folder, 'raw_data')
    #
    # dict_crispr_arrays = import_data.load_crispr_arrays_from_pickle(os.path.join(raw_data_path,
    #                                                                              'dict_cas_gene_sequences.pickle'))
    # nb_repeats = 0
    # nb_spacer_arrays = 0
    # for repeat_name, ls_spacer_arrays in dict_crispr_arrays.items():
    #     if len(ls_spacer_arrays) >= 3:
    #         nb_repeats += 1
    #         nb_spacer_arrays += len(ls_spacer_arrays)
    #         crispr_group = CRISPRGroup(repeat_name, ls_spacer_arrays)
    #
    #         # create_cas_tree(crispr_group, work_path, save_path)
    # print('Nb repeats: ', nb_repeats, 'Nb spacer arrays: ', nb_spacer_arrays)
