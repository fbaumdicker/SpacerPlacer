import json

from model.helpers.misc import create_logger

def reverse_complement(dna_seq):
    """
    :param dna_seq:
    :return:
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join([complement[base] if base in complement else base for base in dna_seq[::-1]])

class CRISPRFinderParser:
    def __init__(self, crisprfinder_result_json_file, min_evidence_level=4, reorient_by_direction_pred=True, logger=None,
                 cluster_by_spacer_overlap=True):
        self.crisprfinder_result_json_file = crisprfinder_result_json_file
        self.min_evidence_level = min_evidence_level
        self.reorient_by_direction_pred = reorient_by_direction_pred
        self.crispr_dictionary = {}
        self.clustered_by_repeat_arrays = {}
        self.clustered_groups = {}
        self.clustered_groups_spacer_ids = {}
        self.cluster_spacer_seq_to_spacer_id = {}
        self.sample_summary = {}
        self.logger = logger if logger is not None else create_logger('CRISPRFinderParser', 1)
        self.cluster_by_spacer_overlap = cluster_by_spacer_overlap

        self._parse_ls_crisprfinder_files()
        self._cluster_array_groups_by_repeat()
        if self.cluster_by_spacer_overlap:
            self._cluster_groups_by_spacer_overlap()
        self._distribute_spacer_ids()
        # self.write_spacer_fasta()

    def write_spacer_fasta(self, output_folder):
        self.logger.info(f'Writing spacer sequences to FASTA files in {output_folder}')
        for group_name, ls_array_tuples in self.clustered_groups_spacer_ids.items():
            with open(f'{output_folder}/{group_name}.fa', 'w') as f:
                for array_name, spacers, direction in ls_array_tuples:
                    f.write(f'>{array_name} ' + ' ; ' + f'{direction}' + ' \n')
                    f.write(f'{", ".join([str(s) for s in spacers])} \n')

    def get_clustered_groups_spacer_ids(self):
        return self.clustered_groups_spacer_ids

    def get_cluster_spacer_seq_to_spacer_id(self):
        return self.cluster_spacer_seq_to_spacer_id

    def get_cluster_spacer_id_to_spacer_seq(self):
        return {v: k for k, v in self.cluster_spacer_seq_to_spacer_id.items()}

    def _distribute_spacer_ids(self):
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


    def _cluster_array_groups_by_repeat(self):
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

    def _cluster_groups_by_spacer_overlap(self):
        self.logger.info('Clustering CRISPR arrays based on spacer overlap')
        for repeat, array_tuples in self.clustered_by_repeat_arrays.items():
            counter = 0
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
                    repeat_clusters[array_name + '_' + str(counter)] = [(array_name, spacers, direction)]
                    counter += 1

            # add repeat clusters to self.clustered_groups
            for cluster_name, repeat_clusters in repeat_clusters.items():
                self.clustered_groups[cluster_name] = repeat_clusters

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

