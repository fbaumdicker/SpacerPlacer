

def reverse_complement(dna_seq):
    """
    :param dna_seq:
    :return:
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([complement[base] if base in complement else base for base in dna_seq[::-1]])


def dna_sequence_clean(dna):
    dna = dna.upper()
    dna = dna.replace(" ", "")
    dna = dna.replace("\n", "")
    dna = dna.replace("\r", "")
    dna = dna.replace("\t", "")
    return dna


def get_acc_num(acc_num_info):
    elements = acc_num_info.split("_")
    return "_".join(elements[:-1])


def remove_non_level_4(dict_crispr_by_repeats):
    dict_crispr_by_repeats_no_non_level_4 = {}
    for key, value in dict_crispr_by_repeats.items():
        new_list = []
        for cr in value:
            if cr.evidence_level == 4:
                new_list.append(cr)
        if new_list:
            dict_crispr_by_repeats_no_non_level_4[key] = new_list
    return dict_crispr_by_repeats_no_non_level_4


def remove_by_cas_distance_threshold(dict_crispr_by_repeats, threshold_distance=10000,
                                     additionally_include_cas_unique=True,
                                     exclude_cas_candidates=True):
    dict_crispr_by_repeats_close_cas = {}
    for key, value in dict_crispr_by_repeats.items():
        new_list = []
        for cr in value:
            cas_unique = cr.has_unique_cas(exclude_cas_candidates=exclude_cas_candidates) \
                if additionally_include_cas_unique else False
            if cr.is_close(threshold_distance=threshold_distance) or cas_unique:
                new_list.append(cr)
        if new_list:
            dict_crispr_by_repeats_close_cas[key] = new_list
    return dict_crispr_by_repeats_close_cas


def check_all_orientation(dict_crispr_by_repeats):
    dict_crispr_by_repeats_all_orientation = {}
    for key, value in dict_crispr_by_repeats.items():
        orientation_values = [cr.orientation for cr in value]
        dict_crispr_by_repeats_all_orientation[key] = orientation_values
    return dict_crispr_by_repeats_all_orientation


def check_rev_complement_in_the_dict(dict_crispr_by_repeats):
    list_repeats_reversed_strand = []
    list_repeats_inconsistent_orientation = []
    list_not_predicted = []
    for key, value in dict_crispr_by_repeats.items():
        for cr in value:
            if cr.orientation == 2:
                list_repeats_reversed_strand.append(key)
                break
    list_repeats_reversed_strand_rev_complement = [reverse_complement(repeat)
                                                   for repeat in list_repeats_reversed_strand]

    list_present_in_dict = []
    for repeat in list_repeats_reversed_strand_rev_complement:
        if repeat in dict_crispr_by_repeats:
            value = dict_crispr_by_repeats[repeat]
            first_cr = value[0]
            if first_cr.orientation == 1:
                list_present_in_dict.append(repeat)
            elif first_cr.orientation == 2:
                list_repeats_inconsistent_orientation.append(repeat)
            else:
                list_not_predicted.append(repeat)
    return list_present_in_dict, list_repeats_inconsistent_orientation, list_not_predicted


def fix_orientation(dict_crispr_by_repeats, dict_spacers_id):
    dict_id_spacers = {value: key for key, value in dict_spacers_id.items()}
    dict_crispr_fixed_orientation = {}
    for key, value in dict_crispr_by_repeats.items():
        orientation = value[0].orientation
        if orientation in [1, 3]:
            if key not in dict_crispr_fixed_orientation:
                dict_crispr_fixed_orientation[key] = value
            else:
                dict_crispr_fixed_orientation[key].extend(value)
        elif orientation == 2:
            key_repeat = reverse_complement(key)
            new_cr_list = []
            for cr in value:
                new_cr = cr
                new_cr.orientation = 1
                new_cr.drconsensus = reverse_complement(cr.drconsensus)
                spacer_indexes = new_cr.spacer_indexes
                spacer_sequences = [dict_id_spacers[spacer_index] for spacer_index in spacer_indexes]
                spacer_sequences = [reverse_complement(spacer) for spacer in spacer_sequences][::-1]
                new_cr.spacer_indexes = [dict_spacers_id[spacer] for spacer in spacer_sequences]
                new_cr_list.append(new_cr)
            if key_repeat in dict_crispr_fixed_orientation:
                dict_crispr_fixed_orientation[key_repeat].extend(new_cr_list)
            else:
                dict_crispr_fixed_orientation[key_repeat] = new_cr_list
    return dict_crispr_fixed_orientation


