import json
import os.path

from sp_model.helpers.misc import create_logger


def set_default(obj):
    if isinstance(obj, set):
        return list(obj)
    raise TypeError


def reverse_complement(dna_seq):
    """
    :param dna_seq:
    :return:
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([complement[base] if base in complement else base for base in dna_seq[::-1]])

def levenshtein_distance(seq1, seq2):
    len_seq1 = len(seq1) + 1
    len_seq2 = len(seq2) + 1

    # Create a matrix
    matrix = [[0] * len_seq2 for _ in range(len_seq1)]

    # Initialize the matrix
    for i in range(len_seq1):
        matrix[i][0] = i
    for j in range(len_seq2):
        matrix[0][j] = j

    # Compute the Levenshtein distance
    for i in range(1, len_seq1):
        for j in range(1, len_seq2):
            if seq1[i - 1] == seq2[j - 1]:
                cost = 0
            else:
                cost = 1
            matrix[i][j] = min(matrix[i - 1][j] + 1,  # Deletion
                               matrix[i][j - 1] + 1,  # Insertion
                               matrix[i - 1][j - 1] + cost)  # Substitution

    return matrix[-1][-1]

class CRISPRFinderParser:
    def __init__(self, crisprfinder_result_json_file, output_folder=None, output_folder_for_too_small_groups=None,
                 min_evidence_level=4,
                 reorient_by_direction_pred=True, logger=None,
                 cluster_by_spacer_overlap=True, filter_more_than_x_in_group=None, max_levenshtein_distance_spacers=0):
        self.crisprfinder_result_json_file = crisprfinder_result_json_file
        self.min_evidence_level = min_evidence_level
        self.reorient_by_direction_pred = reorient_by_direction_pred
        self.max_levenshtein_distance_spacers = max_levenshtein_distance_spacers
        self.crispr_dictionary = {}
        self.clustered_by_repeat_arrays = {}
        self.clustered_groups = {}
        self.clustered_groups_spacer_ids = {}
        self.cluster_spacer_seq_to_spacer_id = {}
        self.spacer_id_to_clustered_spacer_seq_for_each_group = {}
        self.spacer_cluster = None
        self.sample_summary = {}
        self.logger = logger if logger is not None else create_logger('CRISPRFinderParser', 1)
        self.cluster_by_spacer_overlap = cluster_by_spacer_overlap
        self.filter_more_than_x_in_group = filter_more_than_x_in_group
        self.output_folder = output_folder
        self.output_folder_for_too_small_groups = output_folder_for_too_small_groups

        self._parse_ls_crisprfinder_files()
        self._cluster_array_groups_by_repeat()
        if max_levenshtein_distance_spacers > 0:
            self._cluster_spacers_by_levenshtein_distance()
        if self.cluster_by_spacer_overlap:
            self._cluster_groups_by_spacer_overlap()
        self._distribute_spacer_ids()
        if self.output_folder is not None:
            if not os.path.exists(self.output_folder):
                os.makedirs(self.output_folder)
            if not os.path.exists(self.output_folder_for_too_small_groups):
                os.makedirs(self.output_folder_for_too_small_groups)
            self.write_spacer_fasta(self.output_folder, filter_more_than_x_in_group=self.filter_more_than_x_in_group)
            # self.write_spacer_seq_to_spacer_id(self.output_folder)


    def write_spacer_fasta(self, output_folder, filter_more_than_x_in_group=None):
        self.logger.info(f'Writing spacer sequences to FASTA files in {output_folder}')
        for group_name, ls_array_tuples in self.clustered_groups_spacer_ids.items():
            of = output_folder
            if filter_more_than_x_in_group is not None:
                if len(ls_array_tuples) <= filter_more_than_x_in_group:
                    self.logger.info(f'Group {group_name} has less than {filter_more_than_x_in_group} and is saved to {self.output_folder_for_too_small_groups}.')
                    of = self.output_folder_for_too_small_groups
            with open(f'{of}/{group_name}.fa', 'w') as f:
                for array_name, spacers, direction in ls_array_tuples:
                    f.write(f'>{array_name} \n')
                    f.write(f'{", ".join([str(s) for s in spacers])} \n')

    def write_spacer_seq_to_spacer_id(self, output_folder):
        self.logger.info(f'Writing dictionary of spacer sequences -> spacer ids to {output_folder}')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        for group_name, dict_spacer_seq_to_spacer_id in self.cluster_spacer_seq_to_spacer_id.items():
            with open(f'{output_folder}/{group_name}_spacer_seq_to_spacer_id.json', 'w') as f:
                json.dump(dict_spacer_seq_to_spacer_id, f)

    def write_spacer_cluster_to_file(self, output_folder):
        """Saves spacer clusters to file. The first spacer in each cluster is the representative spacer."""
        self.logger.info(f'Writing spacer clusters to {output_folder}')
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        if self.spacer_cluster is not None:
            with open(f'{output_folder}/spacer_clusters.json', 'w') as f:
                json.dump(self.spacer_cluster, f, default=set_default)

    def get_spacer_cluster(self):
        return self.spacer_cluster

    def get_clustered_groups_spacer_ids(self):
        return self.clustered_groups_spacer_ids

    def get_cluster_spacer_seq_to_spacer_id(self):
        return self.cluster_spacer_seq_to_spacer_id

    def get_cluster_spacer_id_to_spacer_seq(self):
        return {v: k for k, v in self.cluster_spacer_seq_to_spacer_id.items()}

    def _parse_ls_crisprfinder_files(self):
        if isinstance(self.crisprfinder_result_json_file, list):
            for file_path in self.crisprfinder_result_json_file:
                self._parse_crisprfinder_file(file_path)
        else:
            self._parse_crisprfinder_file(self.crisprfinder_result_json_file)

    def _parse_crisprfinder_file(self, file_path):
        self.logger.info(f'Parsing CRISPRFinder JSON file: {file_path}')
        with open(file_path, 'r') as f:
            data = json.load(f)
            ls_results = data['Sequences']

            for result in ls_results:
                sequence_id = result['Id']
                sequence_length = result['Length']
                sequence_cas = result['Cas']
                self.logger.info(f'Parsing CRISPRFinder result for sequence {sequence_id} with length {sequence_length}, '
                                 f'skipping CRISPRs with evidence level < {self.min_evidence_level}:')

                sequence_crisprs = result['Crisprs']
                if len(sequence_crisprs) == 0:
                    self.logger.info(f'No CRISPRs found for sequence {sequence_id}')
                    continue
                else:
                    self.logger.info(f'Found {len(sequence_crisprs)} CRISPR(s) for sequence {sequence_id}')

                for crispr in sequence_crisprs:
                    if crispr['Evidence_Level'] < self.min_evidence_level:
                        self.logger.debug(f'Skipping CRISPR with evidence level {crispr["Evidence_Level"]}')
                        continue
                    crispr_name = crispr['Name']
                    crispr_start = crispr['Start']
                    crispr_end = crispr['End']

                    self.logger.info(f'Parsing CRISPR array {crispr_name} with {crispr["Spacers"]} spacers,'
                                     f'consensus repeat: {crispr["DR_Consensus"]}, '
                                     f'direction: {crispr["CRISPRDirection"]},'
                                     f'potential orientation: {crispr["Potential_Orientation"]}.')
                    dr_consensus = crispr['DR_Consensus']
                    crispr_direction = crispr['CRISPRDirection']
                    potential_orientation = crispr['Potential_Orientation']

                    repeats = []
                    spacers = []

                    for val in crispr['Regions']:
                        if val['Type'] == 'DR':
                            repeats.append(val['Sequence'])
                        elif val['Type'] == 'Spacer':
                            spacers.append(val['Sequence'])

                    self.crispr_dictionary[crispr_name] = {'start': crispr_start,
                                                           'end': crispr_end,
                                                           'direction': crispr_direction,
                                                           'potential_orientation': potential_orientation,
                                                           'dr_consensus': dr_consensus,
                                                           'repeats': repeats,
                                                           'spacers': spacers}

    def _cluster_array_groups_by_repeat(self):
        """
        Clusters CRISPR arrays based on equal consensus repeat.
        :return:
        """
        self.logger.info('Clustering CRISPR arrays based on equal consensus repeat')
        for array_name, array_data in self.crispr_dictionary.items():
            direction = array_data['direction']
            dr_consensus = array_data['dr_consensus']
            spacers = array_data['spacers']
            # Not entirely sure how the output of CRISPRCasFinder is structured. For now, I assume that the spacers
            # are in the same orientation as in the input sequence (and thus repeats/arrays/spacers
            # need to be reversed for overlap search).
            if self.reorient_by_direction_pred and direction == '-':
                dr_consensus = reverse_complement(dr_consensus)
                spacers = [reverse_complement(s) for s in spacers[::-1]]
            if dr_consensus not in self.clustered_by_repeat_arrays:
                self.clustered_by_repeat_arrays[dr_consensus] = []
            self.clustered_by_repeat_arrays[dr_consensus].append((array_name, spacers,
                                                                  array_data['direction']))


    def _cluster_spacers_by_levenshtein_distance(self):
        """
        Clusters spacers based on Levenshtein distance (within the groups with same repeat).
        For each spacer cluster a representative is chosen which replaces the spacers from the cluster in all arrays
        (within a repeat group). Also generates a dictionary mapping the selected representative mapped to all spacers
        in the cluster (for each repeat group).
        :return:
        """
        self.logger.info(f'Clustering spacers based on Levenshtein distance with max. distance {self.max_levenshtein_distance_spacers}')
        self.spacer_cluster = {}
        new_repeat_array_tuples = {}
        for repeat, array_tuples in self.clustered_by_repeat_arrays.items():
            self.spacer_cluster[repeat] = {}
            new_repeat_array_tuples[repeat] = []
            for (array_name, spacers, direction) in array_tuples:
                new_spacers = []
                for s in spacers:
                    found_in_cluster = False
                    for representative, cluster in self.spacer_cluster[repeat].items():
                        if any([(levenshtein_distance(s, c) <= self.max_levenshtein_distance_spacers) for c in cluster]):
                            self.spacer_cluster[repeat][representative].add(s)
                            new_spacers.append(representative)
                            found_in_cluster = True
                            break
                    if not found_in_cluster:
                        self.spacer_cluster[repeat][s] = {s}
                        new_spacers.append(s)
                new_repeat_array_tuples[repeat].append((array_name, new_spacers, direction))
        self.clustered_by_repeat_arrays = new_repeat_array_tuples
        return

    def _cluster_groups_by_spacer_overlap(self):
        """
        Clusters CRISPR arrays (within repeat groups) to find subgroups with spacer overlap.
        :return:
        """
        self.logger.info('Clustering CRISPR arrays based on spacer overlap')
        counter = 0
        for repeat, array_tuples in self.clustered_by_repeat_arrays.items():
            if len(array_tuples) == 1:
                self.logger.info(f'Only one array with consensus repeat {repeat} found, skipping clustering.')
                # continue

            repeat_clusters = {}
            for (array_name, spacers, direction) in array_tuples:
                found_in_cluster = False
                for (cluster_name, ls_array_tuples) in repeat_clusters.items():
                    if found_in_cluster:
                        break
                    for cluster_array_name, cluster_spacers, cluster_direction in ls_array_tuples:
                        if len(set(spacers) & set(cluster_spacers)) > 0:
                            repeat_clusters[cluster_name].append((array_name, spacers, direction))
                            found_in_cluster = True
                            break


                if not found_in_cluster:
                    repeat_clusters['g_' + str(counter)] = [(array_name, spacers, direction)]
                    counter += 1

            # add repeat clusters to self.clustered_groups
            for cluster_name, repeat_clusters in repeat_clusters.items():
                self.clustered_groups[cluster_name] = repeat_clusters

    def _distribute_spacer_ids(self):
        """
        Assigns an ID (number) to each unique spacer sequence in each group (IDs are not globally unique across groups).
        (If spacers are clustered by overlap the ID stands for the representative spacer (which in turn stands for all
        spacers in the cluster).) A dictionary is created mapping (representative) spacer sequences to spacer IDs.
        :return:
        """
        self.logger.info('Distributing spacer IDs')
        if self.cluster_by_spacer_overlap:
            groups = self.clustered_groups
        else:
            groups = self.clustered_by_repeat_arrays
        for group_name, ls_array_tuples in groups.items():
            dict_spacer_seq_to_spacer_id = {}
            spacer_counter = 0
            ls_tuples = []
            for array_name, spacers, direction in ls_array_tuples:
                spacers_as_id = []
                for spacer_seq in spacers:
                    if spacer_seq not in dict_spacer_seq_to_spacer_id:
                        spacer_id = spacer_counter
                        dict_spacer_seq_to_spacer_id[spacer_seq] = spacer_id
                        spacer_counter += 1
                    else:
                        spacer_id = dict_spacer_seq_to_spacer_id[spacer_seq]
                    spacers_as_id.append(spacer_id)
                ls_tuples.append((array_name, spacers_as_id, direction))
            self.clustered_groups_spacer_ids[group_name] = ls_tuples
            self.cluster_spacer_seq_to_spacer_id[group_name] = dict_spacer_seq_to_spacer_id

            # TODO implement dictionary spacer_id -> representative spacer seq -> all spacers in cluster
            # This is missing the repeat information, could use repeat as group name (+ additional counter for uniqueness)?
            # dict_spacer_id_to_spacer_clustered_seq = {idy: self.spacer_cluster[repeat][seq] for seq, idy in dict_spacer_seq_to_spacer_id.items()}
            # self.spacer_id_to_clustered_spacer_seq_for_each_group[group_name] = dict_spacer_id_to_spacer_clustered_seq
