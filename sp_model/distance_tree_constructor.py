from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceCalculator, DistanceTreeConstructor
import numpy as np
from scipy.stats import poisson

from Bio.Phylo.Newick import Clade, Tree
import copy

from sp_model import model_tools


MIN_BL = 1e-6
MAX_BL = 20


class LikelihoodDistance(DistanceCalculator):
    """
    Expects arrays to be aligned! But it is not necessary.
    """

    def __init__(self, gain_rate, loss_rate, alpha, provided_lh_fct=None, optimize_separately=True, *args, **kwargs):
        if gain_rate == 0 and loss_rate == 0:
            raise ValueError('Both gain rate and loss rate seem to be 0! No likelihood computation is possible!')
        self.gain_rate = gain_rate
        self.loss_rate = loss_rate
        self.alpha = alpha
        self.provided_lh_fct = provided_lh_fct
        self.big_eps = 1e20
        self.eps = 1e-40
        self.optimize_separately = optimize_separately
        super().__init__(*args, **kwargs)

        self.min_val = MIN_BL
        self.max_val = MAX_BL

    def get_distance(self, msa):
        names = msa[0]
        arrays = msa[1]
        scores = []
        for i, x in enumerate(arrays):
            row = []
            for j, y in enumerate(arrays):
                if j <= i:
                    row.append(self._pairwise(x, y))
            scores.append(row)
        return DistanceMatrix(names, matrix=scores)

    def _pairwise(self, seq1, seq2):
        if seq1 == seq2:
            return 0
        fes_pos, fes = self._find_FES(seq1, seq2)
        ls_gain_counts, ls_keep_counts, ls_deletion_lengths = self._find_gains_deletions(seq1, seq2, fes_pos)
        if self.optimize_separately:
            t_1 = self._optimize_bl_separately(ls_gain_counts[0], ls_keep_counts[0], ls_deletion_lengths[0]).x
            t_2 = self._optimize_bl_separately(ls_gain_counts[1], ls_keep_counts[1], ls_deletion_lengths[1]).x
        else:
            opt_result = self._optimize_bl(ls_gain_counts, ls_keep_counts, ls_deletion_lengths)
            t_1, t_2 = opt_result.x[0], opt_result.x[1]
        t_1_2 = min(t_1 + t_2, self.max_val)
        t_1_2 = max(t_1_2, self.min_val)
        return t_1_2

    def _optimize_bl_separately(self, gain_count, keep_count, deletion_length):
        lh_fct = lambda t: -self._gain_lh([t], [gain_count])[0] \
                            - self._deletion_keep_lh([t], [keep_count], [deletion_length])
        bounds = (self.eps, 1000000)
        return model_tools.minimize_scalar_fct(lh_fct, bounds=bounds, method='bounded')

    def _optimize_bl(self, ls_gain_counts, ls_keep_counts, ls_deletion_lengths):
        lh_fct = lambda t: - sum(self._gain_lh([t[0], t[1]], ls_gain_counts)) \
                           - self._deletion_keep_lh([t[0], t[1]], ls_keep_counts, ls_deletion_lengths)
        bounds = ((self.eps, None), (self.eps, None))
        start_values = np.array([0.1, 0.1])
        return model_tools.minimize_fct(lh_fct, start_values, bounds=bounds, method='L-BFGS-B', tol=1e-10)

    def _gain_lh(self, ls_t, ls_gain_counts):
        lh = [poisson.logpmf(k, t * self.gain_rate / (self.loss_rate * self.alpha + self.eps))
              for t, k in zip(ls_t, ls_gain_counts)]
        return lh

    def _deletion_keep_lh(self, ls_t, ls_keep_counts, ls_deletion_lengths):
        if self.provided_lh_fct is None:
            lh = model_tools.compute_tree_lh(self.loss_rate,
                                             self.alpha,
                                             ls_t,
                                             ls_keep_counts,
                                             ls_deletion_lengths,
                                             log=True)
        else:
            lh = model_tools.compute_tree_lh_for_given_lh_fct(self.loss_rate,
                                                              self.alpha,
                                                              ls_t,
                                                              ls_keep_counts,
                                                              ls_deletion_lengths,
                                                              lambdifyed_lh_fct=self.provided_lh_fct,
                                                              log=True)
        return lh

    def _find_FES(self, seq1, seq2):
        for i in range(max(len(seq1), len(seq2))):
            s_1 = seq1[i]
            s_2 = seq2[i]
            if s_1 in self.skip_letters or s_2 in self.skip_letters:
                continue
            elif s_1 == s_2:
                return i, s_1
        return max(len(seq1), len(seq2)), None

    def _find_gains_deletions(self, seq1, seq2, fes_pos):
        gain_count_1, gain_count_2 = 0, 0
        keep_count_1, keep_count_2 = 0, 0
        ls_deletion_lengths_1, ls_deletion_lengths_2 = [], []
        deletion_length_1 = 0
        deletion_length_2 = 0
        for i in range(max(len(seq1), len(seq2))):
            s_1 = seq1[i]
            s_2 = seq2[i]
            if s_1 in self.skip_letters and s_2 in self.skip_letters:
                continue
            elif s_1 in self.skip_letters:
                if i < fes_pos:
                    gain_count_2 += 1
                else:
                    keep_count_2 += 1
                    deletion_length_1 += 1
                    if deletion_length_2 > 0:
                        ls_deletion_lengths_2.append(deletion_length_2)
                        deletion_length_2 = 0
            elif s_2 in self.skip_letters:
                if i < fes_pos:
                    gain_count_1 += 1
                else:
                    keep_count_1 += 1
                    deletion_length_2 += 1
                    if deletion_length_1 > 0:
                        ls_deletion_lengths_1.append(deletion_length_1)
                        deletion_length_1 = 0
            else:
                if i < fes_pos:
                    # s_1 == s_2: this should normally not happen (should be fes then)
                    # s_1 != s_2: misalignment by mafft -> gain for both of different spacers
                    gain_count_1 += 1
                    gain_count_2 += 1
                else:
                    keep_count_1 += 1
                    keep_count_2 += 1
                    if deletion_length_1 > 0:
                        ls_deletion_lengths_1.append(deletion_length_1)
                        deletion_length_1 = 0
                    if deletion_length_2 > 0:
                        ls_deletion_lengths_2.append(deletion_length_2)
                        deletion_length_2 = 0
        if deletion_length_1 > 0:
            ls_deletion_lengths_1.append(deletion_length_1)
        if deletion_length_2 > 0:
            ls_deletion_lengths_2.append(deletion_length_2)
        return [gain_count_1, gain_count_2], [keep_count_1, keep_count_2], \
            [ls_deletion_lengths_1, ls_deletion_lengths_2]


class FixedDistanceTreeConstructor(DistanceTreeConstructor):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.min_bl = MIN_BL
        self.max_bl = MAX_BL

    def upgma(self, distance_matrix):
        """Construct and return an UPGMA tree.
        Constructs and returns an Unweighted Pair Group Method
        with Arithmetic mean (UPGMA) tree.
        :Parameters:
            distance_matrix : DistanceMatrix
                The distance matrix for tree construction.
        """
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Must provide a DistanceMatrix object.")

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades
        clades = [Clade(None, name) for name in dm.names]
        cardinalities = [1 for _ in dm.names]
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        inner_clade = Clade(None, "Inner")
        while len(dm) > 1:
            min_dist = dm[1, 0]
            # find minimum index
            for i in range(1, len(dm)):
                for j in range(0, i):
                    if min_dist >= dm[i, j]:
                        min_dist = dm[i, j]
                        min_i = i
                        min_j = j

            # create clade
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)

            # assign branch length
            dist = min_dist * 1.0 / 2 - self._height_of(clade1)
            clade1.branch_length = dist

            dist = min_dist * 1.0 / 2 - self._height_of(clade2)
            clade2.branch_length = dist

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (cardinalities[min_i] * dm[min_i, k] + cardinalities[min_j] * dm[min_j, k]) * 1.0 / (
                            cardinalities[min_i] + cardinalities[min_j])

            dm.names[min_j] = "Inner" + str(inner_count)

            cardinalities[min_j] = cardinalities[min_j] + cardinalities[min_i]
            del cardinalities[min_i]

            del dm[min_i]
        inner_clade.branch_length = 0

        return Tree(inner_clade, rooted=True)

    def nj(self, distance_matrix):
        """
        Mostly unchanged from biopython implementation, only added min_bl to avoid negative branch lengths.

        Construct and return a Neighbor Joining tree.

                :Parameters:
                    distance_matrix : DistanceMatrix
                        The distance matrix for tree construction.

                """
        if not isinstance(distance_matrix, DistanceMatrix):
            raise TypeError("Must provide a DistanceMatrix object.")

        # make a copy of the distance matrix to be used
        dm = copy.deepcopy(distance_matrix)
        # init terminal clades
        clades = [Clade(None, name) for name in dm.names]
        # init node distance
        node_dist = [0] * len(dm)
        # init minimum index
        min_i = 0
        min_j = 0
        inner_count = 0
        # special cases for Minimum Alignment Matrices
        if len(dm) == 1:
            root = clades[0]

            return Tree(root, rooted=False)
        elif len(dm) == 2:
            # minimum distance will always be [1,0]
            min_i = 1
            min_j = 0
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            clade1.branch_length = dm[min_i, min_j] / 2.0
            clade2.branch_length = dm[min_i, min_j] - clade1.branch_length
            inner_clade = Clade(None, "Inner")
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            clades[0] = inner_clade
            root = clades[0]

            return Tree(root, rooted=False)
        while len(dm) > 2:
            # calculate nodeDist
            # step 1
            for i in range(0, len(dm)):
                node_dist[i] = 0
                for j in range(0, len(dm)):
                    node_dist[i] += dm[i, j]
                node_dist[i] = node_dist[i] / (len(dm) - 2)
            # find minimum distance pair
            min_dist = dm[1, 0] - node_dist[1] - node_dist[0]
            min_i = 0
            min_j = 1
            # step 2
            for i in range(1, len(dm)):
                for j in range(0, i):
                    temp = dm[i, j] - node_dist[i] - node_dist[j]
                    if min_dist > temp:
                        min_dist = temp
                        min_i = i
                        min_j = j
            # create clade
            # step 3
            clade1 = clades[min_i]
            clade2 = clades[min_j]
            inner_count += 1
            inner_clade = Clade(None, "Inner" + str(inner_count))
            inner_clade.clades.append(clade1)
            inner_clade.clades.append(clade2)
            # assign branch length
            clade1.branch_length = max((dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]) / 2.0, self.min_bl)
            clade2.branch_length = max(dm[min_i, min_j] - clade1.branch_length, self.min_bl)
            # clade1.branch_length = (dm[min_i, min_j] + node_dist[min_i] - node_dist[min_j]) / 2.0
            # clade2.branch_length = dm[min_i, min_j] - clade1.branch_length

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            # step 4
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (dm[min_i, k] + dm[min_j, k] - dm[min_i, min_j]) / 2.0
            dm.names[min_j] = "Inner" + str(inner_count)
            # print('dm', dm)
            del dm[min_i]

        # set the last clade as one of the child of the inner_clade
        root = None
        if clades[0] == inner_clade:
            clades[0].branch_length = 0
            clades[1].branch_length = max(dm[1, 0], self.min_bl)
            clades[0].clades.append(clades[1])
            root = clades[0]
        else:
            clades[0].branch_length = max(dm[1, 0], self.min_bl)
            clades[1].branch_length = 0
            clades[1].clades.append(clades[0])
            root = clades[1]

        return Tree(root, rooted=False)

    def _height_of(self, clade):
        if clade.is_terminal():
            height = 0
        else:
            height = max(c.branch_length + self._height_of(c) for c in clade.clades)
        return height

