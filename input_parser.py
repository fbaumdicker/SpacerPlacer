import os


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

class InputParser:
    def __init__(self, input_folder, restructured_output_path, spacer_number_to_seq_file=None,
                 cluster_spacers=False, check_orientation_for_highest_spacer_overlap=False,
                 split_into_groups_by_spacer_overlap=False):
        self.input_folder = input_folder
        self.restructured_output_path = restructured_output_path
        self.spacer_number_to_seq_file = spacer_number_to_seq_file
        self.spacer_seq_to_number = None
        # self.spacer_number_to_seq = None
        self.file_name_to_spacer_indices = None
        self.spacer_orientations = None

        self.cluster_spacers = cluster_spacers
        self.check_orientation_for_highest_spacer_overlap = check_orientation_for_highest_spacer_overlap
        self.split_into_groups_by_spacer_overlap = split_into_groups_by_spacer_overlap

        self._enumerate_spacers()
        self.clustering_steps()
        self.write_output_to_file()

    def _enumerate_spacers(self):
        spacers = {}
        spacers_ori = {}
        spacer_count = 0
        result = {}

        forward_orientation = os.path.join(self.input_folder, "pos_strand")
        reverse_orientation = os.path.join(self.input_folder, "neg_strand")

        for file_name_index, file_name in enumerate(os.listdir(forward_orientation), 1):
            result[file_name.split('.')[0]] = []
            with open(os.path.join(forward_orientation, file_name), 'r') as file:
                for line in file:
                    if line.startswith(">"):
                        continue
                    spacer = line.strip()
                    if spacer not in spacers:
                        spacers[spacer] = spacer_count
                        spacers_ori[spacer] = "pos_strand"
                        spacer_count += 1
                    result[file_name.split('.')[0]].append(spacers[spacer])

        for file_name_index, file_name in enumerate(os.listdir(reverse_orientation), len(result) + 1):
            result[file_name.split('.')[0]] = []
            with open(os.path.join(reverse_orientation, file_name), 'r') as file:
                lines = file.readlines()
                lines = lines[::-1]
                for line in lines:
                    if line.startswith(">"):
                        continue
                    spacer = line.strip()
                    if self._rev_com(spacer) not in spacers:
                        spacers[self._rev_com(spacer)] = spacer_count
                        spacers_ori[self._rev_com(spacer)] = "neg_strand"
                        spacer_count += 1
                    result[file_name.split('.')[0]].append(spacers[self._rev_com(spacer)])
        self.spacer_orientations = spacers_ori
        self.spacer_seq_to_number = spacers
        # self.spacer_number_to_seq = {v: k for k, v in spacers.items()}
        self.file_name_to_spacer_indices = result
        return

    def write_output_to_file(self):
        with open(self.restructured_output_path, 'w') as file:
            for file_name, spacer_indices in self.file_name_to_spacer_indices.items():
                file.write(f">{file_name}\n")
                file.write(", ".join(f"{index}" for index in spacer_indices))
                file.write("\n")

        if self.spacer_number_to_seq_file is not None:
            with open(self.spacer_number_to_seq_file, 'w') as file:
                for spacer_seq, spacer_number in self.spacer_seq_to_number.items():
                    file.write(f">{spacer_number}, {self.spacer_orientations[spacer_seq]}\n")
                    file.write(f"{spacer_seq}")
                    file.write("\n")

    def clustering_steps(self):
        """
        :return:
        """
        if self.check_orientation_for_highest_spacer_overlap:
            # Accept the orientation as is, as a start. Then check for each array, if reversing would increase overlap.
            # If so, reverse the array. If not, keep as is. Run whole ls of arrays this way one after the other.
            # Could be problematic, if the arrays are split in orientation. But probably most efficient.
            # Can cluster spacers in each step. But probably sufficient to cluster at the end.
            raise NotImplementedError
        if self.cluster_spacers:
            # cluster spacers with levenstein distance. Start with single clusters. Then merge clusters with some
            # distance (1) in mutations iteratively, until no clusters are mergeable.
            new_clusters = self.cluster_spacers_with_levenshtein()
            self.rename_spacers(new_clusters)
        if self.split_into_groups_by_spacer_overlap:
            # Split arrays, if they have no overlap in spacers. Need extra cases for no overlap. Throw warning and then
            # main whole group (w/o split)? Or main only overlapping groups and have additional group without overlap?
            raise NotImplementedError
        return

    def cluster_spacers_with_levenshtein(self):
        spacers = [[spacer] for spacer in self.spacer_seq_to_number.keys()]

        def merge_clusters(cluster1, cluster2):
            for spacer1 in cluster1:
                for spacer2 in cluster2:
                    if levenshtein_distance(spacer1, spacer2) == 1:
                        return cluster1 + cluster2
            return None
        def merge_all_clusters(clusters):
            for i in range(len(clusters)):
                for j in range(i+1, len(clusters)):
                    new_cluster = merge_clusters(clusters[i], clusters[j])
                    if new_cluster is not None:
                        clusters[i] = new_cluster
                        clusters.pop(j)
                        return True
            return False

        while merge_all_clusters(spacers):
            pass

        return spacers

    def rename_spacers(self, clusters):
        # Create new numbering by identifying each cluster with a number
        old_numbering = {spacer: number for number, spacer in self.spacer_seq_to_number.items()}

        new_numbering = {}
        for cluster_number, cluster in enumerate(clusters):
            for spacer in cluster:
                new_numbering[spacer] = cluster_number
        self.file_name_to_spacer_indices = {file_name: [new_numbering[old_numbering[spacer]] for spacer in spacers]
                                            for file_name, spacers in self.file_name_to_spacer_indices.items()}
        self.spacer_seq_to_number = new_numbering


    @staticmethod
    def _rev_com(sequence):
        dict_nuc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join([dict_nuc.get(nuc, "N") for nuc in sequence[::-1]])


