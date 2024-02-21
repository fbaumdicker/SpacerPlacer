import os
import pickle

from additional_data.additional_scripts import iqtree_integration


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


def write_folder_system(path_to_pickled_data, folder_system_path):
    """
    Creates folder system containing
    :param path_to_pickled_data:
    :return:
    """
    with open(path_to_pickled_data, 'rb') as f:
        dict_crispr_groups = pickle.load(f)

    dict_acc_nums_multiple_total = {}
    dict_acc_nums_multiple_per_group = {}
    dict_acc_nums_multiple_different_per_group = {}
    total_number_arrays = 0
    unique_acc_nums = set()
    for group_name, crispr_group in dict_crispr_groups.items():
        subfolder_path = os.path.join(folder_system_path, group_name)
        total_number_arrays += len(crispr_group.get_ls_acc_nums())
        dict_dups = give_count_duplicates_in_list(crispr_group.get_ls_acc_nums())
        if dict_dups:
            dict_acc_nums_multiple_per_group[group_name] = dict_dups
            for acc_num, count in dict_acc_nums_multiple_per_group[group_name].items():
                dict_acc_nums_multiple_total[acc_num] = dict_acc_nums_multiple_total.get(acc_num, 0) + count
        set_acc_nums = set(crispr_group.get_ls_acc_nums())
        unique_acc_nums = unique_acc_nums.union(set_acc_nums)

        if len(set_acc_nums) > 2:
            if not os.path.exists(subfolder_path):
                os.makedirs(subfolder_path)
            with open(os.path.join(subfolder_path, 'acc_nums.txt'), 'w') as f:
                for acc_num in set_acc_nums:
                    f.write(acc_num)
                    f.write('\n')
        else:
            print('Skipping group', group_name, 'because it has less than 3 acc_nums')
        if dict_dups:
            ls = []
            for key in dict_acc_nums_multiple_per_group[group_name].keys():
                ls_arrays = crispr_group.get_ls_arrays_by_acc_num(key)
                print(key)
                for array in ls_arrays:
                    print(array.spacer_array)
                ls.append(all([ls_arrays[0].spacer_array == x.spacer_array for x in ls_arrays]))
            dict_acc_nums_multiple_different_per_group[group_name] = (sum(ls), len(ls))

    print('Total number of arrays: ', total_number_arrays)
    print('Number of unique accession numbers: ', len(unique_acc_nums))
    print('Number of accession numbers with multiple arrays: ', len(dict_acc_nums_multiple_total.keys()))
    print('Maximum number of arrays for one acc num: ', max(dict_acc_nums_multiple_total.values()))
    print('Number of groups with multiple array for one acc num: ', len(dict_acc_nums_multiple_per_group.keys()))
    print('Dict of groups with multiple arrays per acc num: ', dict_acc_nums_multiple_per_group)
    print('Dict of groups with multiple arrays per acc num with all arrays being different '
          '(nb_acc_nums all arrays equal, nb_acc_num with multiple arrays): ',
          dict_acc_nums_multiple_different_per_group)
    print('Total dict of multiple arrays per acc num: ', dict_acc_nums_multiple_total)
    return

def run_iqtree_on_concatenated_core_genomes():
    iqtree_integration.run_iqtree(os.path.join('data', 'concatenated_core_gene_alignments.fa'))
    # tree = import_data.load_single_tree(os.path.join('data', 'concatenated_core_gene_alignments.fa.treefile'))

if __name__ == '__main__':
    path_to_pickled_data = os.path.join('data', '2211_alex_dataset', 'aligned_spacers_data',
                                        'aligned_8_all_reindexed_data_groot',
                                        'dict_aligned_groups.pickle')
    folder_system_path = os.path.join('data', 'v2_folder_system_for_gene_trees', 'aligned_8_all_reindexed_data_groot')
    write_folder_system(path_to_pickled_data, folder_system_path)
