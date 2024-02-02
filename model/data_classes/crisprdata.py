import copy
import numpy as np


class CRISPRGroup:
    def __init__(self, repeat, ls_spacer_arrays, name=None, tree=None):
        self.name = name
        self.repeat = repeat
        ls_name = [spacer_array.name for spacer_array in ls_spacer_arrays]
        self.crispr_dict = {name: spacer_array for name, spacer_array in zip(ls_name, ls_spacer_arrays)}
        self.tree = tree
        # Need to check, if cas types are the same for all arrays in a group.
        self.cas_type = self.determine_group_cas_type()[0]
        self.chromosome_plasmid = self.determine_chromosome_plasmid()

        self.group_species = self.get_group_species()
        self.kingdom = self.get_kingdom()

    def drop_array_by_acc_num(self, acc_num):
        self.crispr_dict = {key: value for key, value in self.crispr_dict.items() if value.acc_num != acc_num}
        return self.crispr_dict

    def get_ls_arrays_by_acc_num(self, acc_num):
        return [a for a in self.crispr_dict.values() if a.acc_num == acc_num]

    def get_ls_acc_nums(self):
        ls_acc_nums = []
        for c_a in self.crispr_dict.values():
            ls_acc_nums.append(c_a.acc_num)
        return ls_acc_nums

    def get_kingdom(self):
        group_kingdom = [a.kingdom for a in self.crispr_dict.values()]
        if len(set(group_kingdom)) > 1:
            return ' + '.join(set(group_kingdom))
        else:
            return group_kingdom[0]

    def get_group_species(self):
        group_species = [a.species_fc for a in self.crispr_dict.values()]
        if len(set(group_species)) > 1:
            return group_species
        else:
            return group_species[0]

    def get_arrays_as_lists(self):
        return {name: crispr_array.get_spacer_array() for name, crispr_array in self.crispr_dict.items()}

    def update_spacer_arrays_by_ls_arrays_as_list(self, ls_array_names, ls_arrays_as_list, remove_other=False):
        for name, array_as_list in zip(ls_array_names, ls_arrays_as_list):
            self.update_spacer_array_by_array_as_list(name, array_as_list)
        if remove_other:
            new_dict_spacer_arrays = {}
            for name in ls_array_names:
                new_dict_spacer_arrays[name] = self.crispr_dict[name]
            self.crispr_dict = new_dict_spacer_arrays
        return self

    def update_spacer_array_by_array_as_list(self, array_name, array_as_list):
        self.crispr_dict[array_name].update_array(array_as_list)
        return

    def determine_group_cas_type(self):
        """
        Computes the majority cas type of group. Will not recognize "similar" types. Returns only one of the max
        count values.
        :return:
        """
        ls_cas_types = list(self.get_cas_type().values())
        if ls_cas_types.count(None) == len(ls_cas_types):
            return None, 0
        nb_types = len(np.unique(ls_cas_types))
        cas_type = max(ls_cas_types, key=ls_cas_types.count)
        return cas_type, nb_types

    def get_cas_type(self):
        dict_cas_type = {}
        for name, spacer_array in self.crispr_dict.items():
            cas_type = spacer_array.get_cas_type()
            dict_cas_type[name] = cas_type
        return dict_cas_type

    def get_cas_gene_info(self):
        dict_cas_info = {}
        for name, spacer_array in self.crispr_dict.items():
            cas_info = spacer_array.all_cas_types()
            dict_cas_info[name] = cas_info
        return dict_cas_info

    def set_tree(self, tree):
        self.tree = tree
        return

    def pop_cas_sequences(self):
        dict_cas_seq = {}
        for name, spacer_array in self.crispr_dict.items():
            cas_seq = spacer_array.pop_cas_sequences()
            dict_cas_seq[name] = cas_seq
        return dict_cas_seq

    def convert_arrays_from_mafft_fmt(self):
        for spacer_array in self.crispr_dict.values():
            spacer_array.convert_arrays_from_mafft_fmt()

    def determine_chromosome_plasmid(self):
        chromosome = 0
        plasmid = 0
        for array in self.crispr_dict.values():
            if array.chromosome_plasmid is None:
                continue
            if 'chromosome' in array.chromosome_plasmid.lower():
                chromosome += 1
            elif 'plasmid' in array.chromosome_plasmid.lower():
                plasmid += 1
            # else:
            #   print('Not chromosome nor plasmid: ', array.name, ' ', array.chromosome_plasmid)
            # print(array.chromosome_plasmid)
        if chromosome > 0 and plasmid > 0:
            label = 'c + p'
        elif chromosome > 0:
            label = 'c'
        elif plasmid > 0:
            label = 'p'
        else:
            label = 'not available'
        return label


class CRISPRArray:
    def __init__(self, name, acc_num, orientation, evidence_level, drconsensus, chromosome_plasmid,
                 kingdom, organism_name, spacer_array, cas_type, cas_genes=None, all_cas_types=None, all_cas_genes=None,
                 species=None, species_fc=None, array_start=None, array_end=None):
        self.acc_num = acc_num
        # We generally expect a dictionary now, probably contains 'direction_prediction', 'strand_prediction',
        # 'strand_confidence'
        self.orientation = orientation
        self.evidence_level = evidence_level
        self.drconsensus = drconsensus
        self.chromosome_plasmid = chromosome_plasmid
        self.kingdom = kingdom
        self.organism_name = organism_name

        self.name = name
        self.spacer_array = [str(s) for s in spacer_array]

        self.cas_type = cas_type
        self.cas_genes = cas_genes

        self.all_cas_types = all_cas_types
        self.all_cas_genes = all_cas_genes

        self.species = species
        # first components of the species tag
        self.species_fc = species_fc
        self.array_start = int(array_start) if array_start is not None else None
        self.array_end = int(array_end) if array_end is not None else None
        self.distance = None

        self._compute_distance()

    def get_spacer_array(self):
        return self.spacer_array

    def has_unique_cas(self, exclude_cas_candidates=True):
        if exclude_cas_candidates:
            return len([t for t in self.all_cas_types if t is not None])
        else:
            return len(self.all_cas_types) == 1

    def is_close(self, threshold_distance=10000):
        if not self.distance:
            return False
        return self.distance <= threshold_distance

    def _compute_distance(self):
        if self.array_start is None or self.array_end is None:
            self.distance = None
            return
        list_cas_intervals = [(cas_gene[1], cas_gene[1] + cas_gene[2]) for cas_gene in self.cas_genes]
        flat_list_cas_intervals = [item for sublist in list_cas_intervals for item in sublist]
        dictances_to_cas_start = [abs(x - self.array_start) for x in flat_list_cas_intervals]
        dictances_to_cas_end = [abs(x - self.array_end) for x in flat_list_cas_intervals]
        combined_distances = dictances_to_cas_start + dictances_to_cas_end
        if not combined_distances:
            self.distance = None
        else:
            self.distance = min(combined_distances)

    def __str__(self):
        return f"<<Name: {self.name}, Acc_number: {self.acc_num}, Orientation: {self.orientation}, " \
               f"Evidence_level: {self.evidence_level}" \
               f", Consensus repeat: {self.drconsensus}" \
               f", Chromosome/Plasmid: {self.chromosome_plasmid}, Kingdom: {self.kingdom}, Name: {self.organism_name}" \
               f", Spacer array: {self.spacer_array}" \
               f", Cas type: {self.cas_type}, All Cas types: {self.all_cas_types}" \
               f", Cas gene sequences {self.cas_genes}>>"

    def get_cas_type(self):
        return self.cas_type

    def update_array(self, array_as_list):
        self.spacer_array = array_as_list

    def all_cas_types(self):
        return self.all_cas_types

    def pop_cas_sequences(self):
        cas_gene_sequences = copy.deepcopy(self.cas_genes)
        all_cas_genes = copy.deepcopy(self.all_cas_genes)
        self.cas_genes = None
        self.all_cas_genes = None
        return cas_gene_sequences, all_cas_genes

    def convert_arrays_from_mafft_fmt(self):
        self.spacer_array = ['-' if s == '--' else s for s in self.spacer_array]
        return

# Data that takes long for no reason:
# GGCTCATCCCCGCTGGCGCGGGGAGCAC
