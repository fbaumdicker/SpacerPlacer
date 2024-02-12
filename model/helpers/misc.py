import itertools

import numpy as np
from timeit import default_timer as timer
import logging
import sys

from Bio import Phylo

from types import ModuleType, FunctionType
from gc import get_referents

BLACKLIST = type, ModuleType, FunctionType


def getsize(obj):
    """sum size of object & members."""
    if isinstance(obj, BLACKLIST):
        raise TypeError('getsize() does not take argument of type: ' + str(type(obj)))
    seen_ids = set()
    size = 0
    objects = [obj]
    while objects:
        need_referents = []
        for obj in objects:
            if not isinstance(obj, BLACKLIST) and id(obj) not in seen_ids:
                seen_ids.add(id(obj))
                size += sys.getsizeof(obj)
                need_referents.append(obj)
        objects = get_referents(*need_referents)
    return size


def sizeof_fmt(num, suffix='B'):
    ''' by Fred Cirera,  https://stackoverflow.com/a/1094933/1870254, modified'''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f %s%s" % (num, 'Yi', suffix)


def create_logger(name, level, outfile=None):
    logger = logging.getLogger(name)
    logger.setLevel(level)

    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(level)
    # formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    ch.setFormatter(formatter)

    logger.addHandler(ch)

    if outfile:
        fh = logging.FileHandler(outfile)

        fh.setLevel(level)
        fh.setFormatter(formatter)

        logger.addHandler(fh)
    return logger


class RunTimer:
    def __init__(self):
        self.start_time = timer()
        self.checkpoints = [self.start_time]

    def time_from_start(self):
        self.checkpoints.append(timer())
        return self.checkpoints[-1] - self.start_time

    def time_between_checkpoints(self, i, j):
        interval = self.checkpoints[j] - self.checkpoints[i]
        return interval

    def time_from_last_checkpoint(self):
        self.checkpoints.append(timer())
        interval = self.checkpoints[-1] - self.checkpoints[-2]
        return interval

    def __str__(self):
        time_from_start = [self.time_between_checkpoints(i, 0) for i in range(len(self.checkpoints))]
        return f'Time since start of Timer {self.checkpoints[-1] - self.start_time};  ' \
               f'Checkpoint times (wrt. start): {time_from_start}'


def remove_completely_same_arrays(ls_arrays, ls_names):
    ls_cleaned_arrays = []
    set_idx_to_remove = set([])
    ls_cleaned_names = []

    dict_combined_names = {}
    for (i, x), (j, y) in itertools.combinations(enumerate(ls_arrays), 2):
        if x == y:
            if j not in set_idx_to_remove:
                dict_combined_names[ls_names[i]] = dict_combined_names.get(ls_names[i], []) + [ls_names[j]]
            set_idx_to_remove.add(j)
    for i, (x, name) in enumerate(zip(ls_arrays, ls_names)):
        if i not in set_idx_to_remove:
            ls_cleaned_arrays.append(x)
            ls_cleaned_names.append(name)

    return ls_cleaned_arrays, ls_cleaned_names, dict_combined_names


class GraphContainsCycleException(Exception):
    pass


def topological_sort(s, adj_matrix, label_dict):
    sorted_spacers = []
    while s:
        n = s.pop()
        sorted_spacers.append(n)
        for j in range(adj_matrix.shape[1]):
            if adj_matrix[n, j] == 1:
                adj_matrix[n, j] = 0
                if sum(adj_matrix[:, j]) == 0:
                    s += [j]
    if np.all(adj_matrix == 0):
        return [label_dict[i] for i in sorted_spacers]
    else:
        raise GraphContainsCycleException(f'Spacer arrays contain (at least one) cycle!')


def estimate_gr_based_on_estimated_lr(ls_avg_array_len, ls_loss_rate, deletion_model, ls_alpha=None):
    ls_gain_rate = []
    for (i, length), lr in zip(enumerate(ls_avg_array_len), ls_loss_rate):
        if deletion_model == 'independent':
            ls_gain_rate.append(lr * length)
        elif deletion_model == 'block':
            ls_gain_rate.append(lr * length * ls_alpha[i])
        else:
            raise NotImplementedError('Estimating gr is not implemented for this model!')
    return ls_gain_rate


def write_tree_to_file(tree, path, fmt='newick'):
    with open(path, 'w') as file:
        Phylo.write(tree, file, format=fmt)
    return


def keep_len_constant(ls_gain_rate, ls_avg_array_length, ls_deletion_model, ls_alpha=None):
    ls_loss_rate = []
    for (i, gr), avg_len, del_m in zip(enumerate(ls_gain_rate), ls_avg_array_length, ls_deletion_model):
        if del_m == 'independent':
            ls_loss_rate.append(gr / avg_len)
        elif del_m == 'block':
            ls_loss_rate.append(gr / (avg_len * ls_alpha[i]))
        else:
            raise NotImplementedError('Keeping lr constant is not implemented for this model!')
    return ls_loss_rate


def uniqify_ls_order_preserving(seq, idfun=None):
    # order preserving
    if idfun is None:
        def idfun(x): return x
    seen = {}
    result = []
    for item in seq:
        marker = idfun(item)
        if marker in seen: continue
        seen[marker] = 1
        result.append(item)
    return result
