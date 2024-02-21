import argparse
import os


class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


class InputParser:
    def __init__(self, input_folder, restructured_input_file, spacer_number_to_seq_file=None):
        self.input_folder = input_folder
        self.restructured_input_file = restructured_input_file
        self.spacer_number_to_seq_file = spacer_number_to_seq_file

        self.enumerate_spacers_and_create_file()

    def enumerate_spacers_and_create_file(self):
        spacers = {}
        spacers_ori = {}
        spacer_count = 1
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

        with open(self.restructured_input_file, 'w') as file:
            for file_name, spacer_indices in result.items():
                file.write(f">{file_name}\n")
                file.write(", ".join(f"{index}" for index in spacer_indices))
                file.write("\n")

        if self.spacer_number_to_seq_file is not None:
            with open(self.spacer_number_to_seq_file, 'w') as file:
                for spacer_seq, spacer_number in spacers.items():
                    file.write(f">{spacer_number}, {spacers_ori[spacer_seq]}\n")
                    file.write(f"{spacer_seq}")
                    file.write("\n")

    @staticmethod
    def _rev_com(sequence):
        dict_nuc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
        return ''.join([dict_nuc.get(nuc, "N") for nuc in sequence[::-1]])


