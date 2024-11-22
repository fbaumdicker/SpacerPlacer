import json

from additional_data.additional_scripts.model.helpers.raw_data_tools import reverse_complement
from model.helpers.misc import create_logger


class CRISPRFinderParser:
    def __init__(self, crisprfinder_result_json_file, min_evidence_level=4, reorient_by_direction_pred=True, logger=None):
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

        self._parse_crisprfinder_file()
        self._cluster_array_groups_by_repeat()
        self._cluster_groups_by_spacer_overlap()
        self._distribute_spacer_ids()

    def get_clustered_groups_spacer_ids(self):
        return self.clustered_groups_spacer_ids

    def get_cluster_spacer_seq_to_spacer_id(self):
        return self.cluster_spacer_seq_to_spacer_id

    def get_cluster_spacer_id_to_spacer_seq(self):
        return {v: k for k, v in self.cluster_spacer_seq_to_spacer_id.items()}

    def _distribute_spacer_ids(self):
        for group_name, ls_array_tuples in self.clustered_groups.items():
            self.cluster_spacer_seq_to_spacer_id[group_name] = {}
            spacer_counter = 0
            for array_name, spacers, direction in ls_array_tuples:
                spacers_as_id = []
                for spacer_seq in spacers:
                    if spacer_seq not in self.cluster_spacer_seq_to_spacer_id[group_name]:
                        spacer_id = spacer_counter
                        self.cluster_spacer_seq_to_spacer_id[spacer_seq] = spacer_id
                        spacer_counter += 1
                    else:
                        spacer_id = self.cluster_spacer_seq_to_spacer_id[spacer_seq]
                    spacers_as_id.append(spacer_id)
                self.clustered_groups_spacer_ids[group_name] = (array_name, spacers_as_id, direction)




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
        for repeat, array_tuples in self.clustered_by_repeat_arrays.items():
            counter = 0
            if len(array_tuples) == 1:
                self.logger.info(f'Only one array with consensus repeat {repeat} found, skipping clustering.')
                continue

            repeat_clusters = {}
            for (array_name, spacers, direction) in array_tuples:
                found_in_cluster = False
                for (cluster_name, ls_array_tuples) in repeat_clusters.items():
                    for cluster_array_name, cluster_spacers, cluster_direction in ls_array_tuples:
                        if len(set(spacers) & set(cluster_spacers)) > 0:
                            repeat_clusters[cluster_name].append((array_name, spacers, direction))
                            found_in_cluster = True

                if not found_in_cluster:
                    repeat_clusters[array_name + '_' + str(counter)] = [(array_name, spacers, direction)]
                    counter += 1

            # add repeat clusters to self.clustered_groups
            for cluster_name, repeat_clusters in repeat_clusters.items():
                self.clustered_groups[cluster_name] = repeat_clusters


    def _parse_crisprfinder_file(self):
        self.logger.info(f'Parsing CRISPRFinder JSON file: {self.crisprfinder_result_json_file}')
        with open(self.crisprfinder_result_json_file, 'r') as f:
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

import os
path_to_crisprfinder_folder = '../../data/2410_amanda_usa_hotspring_data/amanda_nanopore_ccf_results'
ls_assemblies = ['Result_ms50_final_1728918132', 'Result_ms55_final_1728978061', 'Result_ms60_final_1729003293',
            'Result_ms65_final_1729061672']
output_path = '0_result_folder'
for file_path in ls_assemblies:
    results_json_path = os.path.join(path_to_crisprfinder_folder, file_path, 'result.json')
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    parser = CRISPRFinderParser(results_json_path)
    clustered_groups = parser.get_clustered_groups()
    for group_name, ls_array_tuples in clustered_groups.items():
        print(f'Group {group_name}:')
        for array_name, spacers, direction in ls_array_tuples:
            print(f'{array_name}: {spacers} ({direction})')
        print('\n')