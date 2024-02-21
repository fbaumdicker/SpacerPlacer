import numpy as np
from Bio.Phylo.Newick import Clade, Tree


class AdvancedTree(Tree):
    """
    Implements some basic functions for trees. If a tree is not rooted, it will be rooted at the midpoint.
    """
    def __init__(self, root, rooted, model_name=None, *args, **kwargs):
        root = root.root if isinstance(root, Tree) else root
        super(AdvancedTree, self).__init__(root=root,
                                           rooted=rooted,
                                           *args,
                                           **kwargs)
        self.parents = self.determine_parents()
        self.model_name = model_name if model_name is not None else 'advanced_tree'
        if not rooted:
            self.root_at_midpoint()
            self.rooted = True

    def determine_parents(self):
        """
        Sets the up (i.e. parent node) parameter for each clade in the tree. {Returns a dict with clade: parent}
        :return:
        """
        parents = {self.root: self.root}
        self.root.up = None
        for clade in self.root.find_clades(order='level'):
            for child in clade:
                parents[child] = clade
                child.up = clade
        return parents

    def total_branch_length(self):
        return sum(node.branch_length for node in self.root.find_clades())

    def min_max_tree_bl(self):
        min_bl = np.inf
        max_bl = -np.inf
        for clade in self.root.find_clades():
            if clade.up is None:
                continue
            min_bl = min(clade.branch_length, min_bl)
            max_bl = max(clade.branch_length, max_bl)
        return min_bl, max_bl

    def nb_leafs(self):
        return len(list(self.root.get_terminals()))

    def distance_to_leafs(self, return_dict_distances=False):
        """
        Returns maximum distance to leaf (i.e. distance for all leafs if ultrametric tree).
        If return_ls_distances returns dictionary {name of leaf: distance}.
        :return:
        """
        leafs = self.root.get_terminals()
        if leafs:
            if return_dict_distances:
                return {leaf.name: self.distance(leaf) for leaf in leafs}
            else:
                return max(self.distance(leaf) for leaf in leafs)
        else:
            return 0

    def is_ultrametric(self):
        """
        Checks if tree is ultrametric.
        :return:
        """
        return len(set(self.distance_to_leafs(return_dict_distances=True).values())) == 1
