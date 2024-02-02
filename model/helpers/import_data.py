import os
import pickle
import Bio.Phylo as Phylo
import io


def unpickle(f):
    """
    Reads complete pkl file from given file handle.
    :param f:
    :return:
    """
    ls = []
    while True:
        try:
            ls.append(pickle.load(f))
        except EOFError:
            break
    return ls


def load_crispr_arrays_from_pickle(path):
    with open(path, "rb") as f:
        dict_crispr_by_repeats = pickle.load(f)
    return dict_crispr_by_repeats


def load_real_data_tree(path, ls_file_names):
    ls_trees = []
    for file_name in ls_file_names:
        ls_trees += Phylo.parse(os.path.join(path, file_name), format='newick')
    return ls_trees


def load_single_tree(path):
    return next(Phylo.parse(path, format='newick'))


def load_single_tree_from_string(string):
    file_handle = io.StringIO()
    file_handle.write(string)
    file_handle.seek(0)
    return next(Phylo.parse(file_handle, 'newick'))


def read_tomida_data(path, ls_files):
    ls_metadata = []
    ls_arrays = []
    for file in ls_files:
        print(os.path.join(path, file))
        with open(os.path.join(path, file)) as f:
            lines = f.readlines()
            head = lines[0][1:-2]
            head = [h.strip() for h in head.split(',')]
            for line in lines[1:]:
                line = line.split(':')
                metadata = line[0]
                metadata = metadata.strip('()').split(',')
                metadata = [m.strip() for m in metadata]
                arrays = line[1]
                arrays = arrays.split(',')
                arrays = [a.strip() for a in arrays]
                ls_metadata.append(metadata)
                ls_arrays.append(arrays)
    return {'head': head, 'metadata': ls_metadata, 'arrays': ls_arrays}


def read_omer_data(path, ls_files):
    ls_metadata = []
    ls_arrays = []
    for file in ls_files:
        print(os.path.join(path, file))
        with open(os.path.join(path, file)) as f:
            lines = (line.rstrip() for line in f.readlines())
            lines = [line for line in lines if line]
            head = lines[0].split(')')[0]
            head = head.strip('(')
            head = [h.strip().strip("'") for h in head.split(',')]
            # print(head)
            for line in lines[1:]:
                line = line.split(':')
                metadata = line[0]
                arrays = line[1]
                arrays = arrays.split(',')
                arrays = [a.strip() for a in arrays]
                ls_metadata.append(metadata)
                ls_arrays.append(arrays)
    return {'head': head, 'metadata': ls_metadata, 'arrays': ls_arrays}


def read_aligned_fasta(file_path):
    ls_arrays = []
    ls_names = []
    with open(file_path, 'r') as f:
        output = f.readlines()
        for o in output:
            if o[0] == '>':
                ls_names.append(o.strip('>\n'))
                converted_array = []
                ls_arrays.append(converted_array)
            else:
                array = o.strip('\n')
                ls_array = array.split(' ')

                for pos in ls_array:
                    if pos == '--':
                        converted_array.append('-')
                    elif pos == '':
                        continue
                    else:
                        converted_array.append(pos)
    return {'metadata': ls_names, 'arrays': ls_arrays}
