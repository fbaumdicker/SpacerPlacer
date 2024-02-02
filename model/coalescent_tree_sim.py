import numpy as np
import pandas as pd
from Bio.Phylo.Newick import Clade


def rnd_coalescent_tree_distance_matrix_string(nb_leafs):
    tree = rnd_coalescent_tree(nb_leafs, rnd_bl=False)
    # tree.ladderize()
    d_m = calculate_distance_matrix(tree)
    pd_dm = pd.DataFrame(d_m)
    return pd_dm.to_string(header=False, index=False)


def coalescent_branch_lengths(lineages, rnd=False):
    if rnd:
        return [np.random.exponential(2 / (k * (k - 1))) for k in range(2, lineages + 1)]
    return [2 / (k * (k - 1)) for k in range(2, lineages + 1)]


def calculate_distance_matrix(clade):
    d_matrix = np.zeros((len(clade.get_terminals()), len(clade.get_terminals())))
    for i, leaf in enumerate(clade.get_terminals()):
        all_terminals = clade.get_terminals()
        for j, other in enumerate(all_terminals):
            d_matrix[i, j] = clade.distance(leaf, other)
    return d_matrix


def rnd_coalescent_tree(nb_leafs, rnd_bl=False, seed=None, tree_name='t'):
    """
    Generates a random coalescent tree.
    :param rnd_bl:
    :param nb_leafs:
    :param seed:
    :param tree_name:
    :return:
    """
    if seed is not None:
        np.random.seed(seed)
    ls_bl = coalescent_branch_lengths(nb_leafs, rnd=rnd_bl)[::-1]
    name_idx = 0
    ls_clades = [Clade(branch_length=0, name=tree_name + str(i)) for i in range(nb_leafs)]
    # I think this is awful and inefficient
    for bl in ls_bl:
        for clade in ls_clades:
            clade.branch_length = clade.branch_length + bl
        children = []
        for _ in range(2):
            idx = np.random.choice(len(ls_clades))
            child = ls_clades.pop(idx)
            children.append(child)

        new_clade = Clade(branch_length=0, clades=children)
        ls_clades.append(new_clade)
    # sort tree by number of terminal nodes (i.e. single leafs on top)
    ls_clades[0].ladderize()
    # naming, what convention should I use??
    for clade in ls_clades[0].get_terminals():
        clade.name = tree_name + str(name_idx)
        name_idx += 1
    for clade in list(ls_clades[0].get_nonterminals(order='level'))[::-1]:
        clade.name = tree_name + str(name_idx)
        name_idx += 1
    ls_clades[0].name = tree_name
    return ls_clades[0]


def mult_rnd_coalescent_tree(nb_trees, nb_leafs, rnd_bl=False, seed=None, tree_name='t', unique=False):
    """
    Returns a list of random coalescent trees. Note: The list can contain the same (or isomorphic, with random branch
    length homeomorphic(?)) trees multiple times.
    :param rnd_bl:
    :param unique:
    :param nb_trees:
    :param nb_leafs:
    :param seed:
    :param tree_name:
    :return:
    """
    if seed is not None:
        np.random.seed(seed)
    ls = []
    while len(ls) < nb_trees:
        if unique:
            for i, tree in enumerate(ls[:-1]):
                if compare(ls[-1], tree):
                    del ls[i]
                    # print(ls)
        ls.append(rnd_coalescent_tree(nb_leafs, rnd_bl=rnd_bl, tree_name=tree_name))
    return ls


def compare(tree1, tree2):
    term_names1 = [term.branch_length for term in tree1.get_terminals()]
    term_names2 = [term.branch_length for term in tree2.get_terminals()]
    # false if terminals are not the same
    if term_names1 != term_names2:
        return False
    else:
        return True
