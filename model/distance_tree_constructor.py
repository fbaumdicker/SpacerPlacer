import pandas as pd
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceCalculator, DistanceTreeConstructor
import numpy as np
from scipy.stats import poisson

from Bio.Phylo.Newick import Clade, Tree
import copy

from model import model_tools

# height excursion explodes for too many arrays
MIN_BL = 1e-6  # 0.000001  # or, if smaller, 0.1/alignment_length in iqtree
MAX_BL = 20  # 2*iqtree value


class LikelihoodDistance(DistanceCalculator):
    """
    Expects arrays to be aligned! But it is not necessary.
    """

    def __init__(self, gain_rate, loss_rate, alpha, provided_lh_fct=None, *args, **kwargs):
        if gain_rate == 0 and loss_rate == 0:
            raise ValueError('Both gain rate and loss rate seem to be 0! No likelihood computation is possible!')
        self.gain_rate = gain_rate
        self.loss_rate = loss_rate
        self.alpha = alpha
        self.provided_lh_fct = provided_lh_fct
        self.big_eps = 1e20
        self.eps = 1e-40
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
        fes_pos, fes = self._find_FES(seq1, seq2)
        ls_gain_counts, ls_keep_counts, ls_deletion_lengths = self._find_gains_deletions(seq1, seq2, fes_pos)
        opt_result = self._optimize_bl(ls_gain_counts, ls_keep_counts, ls_deletion_lengths)
        t_1, t_2 = opt_result.x[0], opt_result.x[1]
        t_1_2 = min(t_1 + t_2, self.max_val)
        t_1_2 = max(t_1_2, self.min_val)
        return t_1_2

    def _optimize_bl(self, ls_gain_counts, ls_keep_counts, ls_deletion_lengths):
        lh_fct = lambda t: - sum(np.log(self._gain_lh([t[0], t[1]], ls_gain_counts))) \
                           - self._deletion_keep_lh([t[0], t[1]], ls_keep_counts, ls_deletion_lengths)
        bounds = ((0.0, None), (0.0, None))
        start_values = np.array([0.1, 0.1])
        return model_tools.minimize_fct(lh_fct, start_values, bounds=bounds, method='Nelder-Mead')

    def _gain_lh(self, ls_t, ls_gain_counts):
        lh = [poisson.pmf(k, t * self.gain_rate / (self.loss_rate * self.alpha + self.eps))
              for t, k in zip(ls_t, ls_gain_counts)]
        # print(lh)
        return [max(self.eps, val) for val in lh]

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
        # dm_counts = copy.deepcopy(distance_matrix)
        # for i in range(1, len(dm_counts)):
        #     for j in range(0, i):
        #         dm_counts[i, j] = 1
        # print(dm_counts)
        # init terminal clades
        clades = [Clade(None, name) for name in dm.names]
        cardinalities = [1 for _ in dm.names]
        # print(cardinalities)
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
            # dist = min(min_dist * 1.0 / 2 - self._height_of(clade1), self.max_bl)
            # dist = max(dist, self.min_bl)
            dist = min_dist * 1.0 / 2 - self._height_of(clade1)
            clade1.branch_length = dist

            # dist = min(min_dist * 1.0 / 2 - self._height_of(clade2), inner_count*self.max_bl)
            # dist = max(dist, self.min_bl)
            dist = min_dist * 1.0 / 2 - self._height_of(clade2)
            clade2.branch_length = dist

            # print('bl', clade1.branch_length, clade2.branch_length)
            # print('heights', self._height_of(clade1), self._height_of(clade2))
            # print('isterimal', clade1.is_terminal(), clade2.is_terminal())
            # if clade1.is_terminal():
            #     clade1.branch_length += self._height_of(clade2)
            # if clade2.is_terminal():
            #     clade2.branch_length += self._height_of(clade1)
            # print('bl after', clade1.branch_length, clade2.branch_length)

            # update node list
            clades[min_j] = inner_clade
            del clades[min_i]

            # print(clades)

            # rebuild distance matrix,
            # set the distances of new node at the index of min_j
            for k in range(0, len(dm)):
                if k != min_i and k != min_j:
                    dm[min_j, k] = (cardinalities[min_i] * dm[min_i, k] + cardinalities[min_j] * dm[min_j, k]) * 1.0 / (
                            cardinalities[min_i] + cardinalities[min_j])

            dm.names[min_j] = "Inner" + str(inner_count)

            cardinalities[min_j] = cardinalities[min_j] + cardinalities[min_i]
            del cardinalities[min_i]
            # print(cardinalities)

            del dm[min_i]
        inner_clade.branch_length = 0

        # for clade in inner_clade.find_clades():
        #     c1 = clade.clades[0]
        #     c2 = clade.clades[1]
        #
        #     if c1.is_terminal():
        #         c1.branch_length += self._height_of(c2)
        #     if c2.is_terminal():
        #         c2.branch_length += self._height_of(c1)
        return Tree(inner_clade)

    def _height_of(self, clade):
        if clade.is_terminal():
            height = 0
        else:
            height = max(c.branch_length + self._height_of(c) for c in clade.clades)
        return height


class NWMSA(DistanceCalculator):
    def __init__(self, match=1000, gap=-.1, mismatch=-10000, skip_letters=None):
        if skip_letters is None:
            skip_letters = ['-']
        super().__init__(skip_letters=skip_letters)
        self.match = match
        self.gap = gap
        self.mismatch = mismatch

    def _pairwise(self, seq1, seq2):
        return 1 - self.nw_score(seq1, seq2) / (self.match * max(len([s for s in seq1 if s not in self.skip_letters]),
                                                                 len([s for s in seq2 if s not in self.skip_letters])))

    def get_distance(self, msa):
        names = msa[0]
        arrays = msa[1]
        scores = []
        for i, x in enumerate(arrays):
            row = []
            for j, y in enumerate(arrays):
                if j <= i:
                    # print(self._pairwise(x, y))
                    # print(x)
                    # print(y)
                    row.append(self._pairwise(x, y))
            scores.append(row)
        return DistanceMatrix(names, matrix=scores)

    def nw_matrix(self, x, y):
        m = np.zeros((len(x) + 1, len(y) + 1))
        m[:, 0] = [a * self.gap for a in range(len(x) + 1)]
        m[0, :] = [a * self.gap for a in range(len(y) + 1)]
        tmp = np.zeros(3)
        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                tmp[0] = m[i - 1, j - 1]
                if x[i - 1] in self.skip_letters and y[j - 1] in self.skip_letters:
                    # tmp[0] += self.gap
                    continue
                else:
                    if x[i - 1] == y[j - 1]:
                        tmp[0] += self.match
                    else:
                        tmp[0] += self.mismatch
                tmp[1] = m[i - 1, j] + self.gap
                tmp[2] = m[i, j - 1] + self.gap
                m[i, j] = max(tmp)
        # print(m)
        return m

    def nw_score(self, x, y):
        return self.nw_matrix(x, y)[-1, -1]

    def get_nw_scores(self, msa):
        names = msa[0]
        arrays = msa[1]
        scores = []
        for i, x in enumerate(arrays):
            row = []
            for j, y in enumerate(arrays):
                if j <= i:
                    row.append(self.nw_score(x, y))
            scores.append(row)
        return DistanceMatrix(names, matrix=scores)

    def nw_alignment(self, x, y):
        aligned_x = []
        aligned_y = []
        m = self.nw_matrix(x, y)
        i, j = len(x), len(y)
        while i > 0 and j > 0:
            s = m[i - 1, j - 1]
            if x[i - 1] in self.skip_letters and y[j - 1] in self.skip_letters:
                s += 0
            elif x[i - 1] == y[j - 1]:
                s += self.match
            else:
                s += self.mismatch
            if i > 0 and j > 0 and m[i, j] == s:
                aligned_x = [x[i - 1]] + aligned_x
                aligned_y = [y[j - 1]] + aligned_y
                i -= 1
                j -= 1
            elif i > 0 and m[i, j] == m[i - 1, j] + self.gap:
                aligned_x = [x[i - 1]] + aligned_x
                aligned_y = ['-'] + aligned_y
                i -= 1
            else:
                aligned_x = ['-'] + aligned_x
                aligned_y = [y[j - 1]] + aligned_y
                j -= 1
        while i > 0:
            aligned_x = [x[i - 1]] + aligned_x
            aligned_y = ['-'] + aligned_y
            i -= 1
        while j > 0:
            aligned_x = ['-'] + aligned_x
            aligned_y = [y[j - 1]] + aligned_y
            j -= 1
        score = m[-1, -1]
        # print('x', x)
        # print('y', y)
        # print('a_x', aligned_x)
        # print('a_y', aligned_y)
        return aligned_x, aligned_y, score

    # def full_alignment(self, tree=None, method=''):
    #     """
    #
    #     :return:
    #     """
    #     if tree is None:
    #         # lengths = []
    #         # for array in self.arrays:
    #         #     lengths.append(len([a for a in array if a not in self.skip_letters]))
    #         # argmax, _ = max(enumerate(lengths), key=lambda x: x[1])
    #         a_x = self.arrays[0]
    #         for y in self.arrays[1:]:
    #             a_x, a_y, _ = self.nw_alignment(a_x, y)
    #             # print(a_x, a_y)
    #             for i, val in enumerate(a_x):
    #                 if val in self.skip_letters:
    #                     a_x[i] = a_y[i]
    #             a_x = [a for a in a_x if a not in self.skip_letters]
    #     else:
    #         for c in tree.find_clades(order='postorder'):
    #             if c.is_terminal():
    #                 c.aligned_array = self.arrays[self.names.index(c.name)]
    #                 c.dup = return_duplicates(c.aligned_array)
    #                 continue
    #             clades = copy.deepcopy(c.clades)
    #             # print(c.clades)
    #             c1 = clades.pop(-1)
    #             a_x = c1.aligned_array
    #             # print('joint node', c.name)
    #             a_y = c1.aligned_array
    #             while len(clades) > 0:
    #                 c2 = clades.pop(-1)
    #                 # print(c1.name, 'c1.aligned arrays', c1.aligned_array)
    #                 # print(c2.name, 'second aligned', c2.aligned_array)
    #                 # print('intersection', set(c1.aligned_array).intersection(set(c2.aligned_array)))
    #                 a_x, a_y, _ = self.nw_alignment(a_x, c2.aligned_array)
    #                 # print('ax', a_x)
    #                 # print('ay', a_y)
    #                 for i, val in enumerate(a_x):
    #                     if val in self.skip_letters:
    #                         a_x[i] = a_y[i]
    #                 # print('ax_full', a_x)
    #
    #                 c.dup = return_duplicates(a_x)
    #                 # print('c.dup', c.dup)
    #                 # print('c2.dup', c2.dup)
    #                 if c2.is_terminal():
    #                     continue
    #                 else:
    #                     diff = len(c.dup) - len(c2.dup)
    #                     if diff > 0:
    #                         a_x_r = self.reversed_alignment(c2)
    #                         a_x_r, a_y_r, _ = self.nw_alignment(c1.aligned_array, a_x_r)
    #
    #                         for i, val in enumerate(a_x_r):
    #                             if val in self.skip_letters:
    #                                 a_x_r[i] = a_y_r[i]
    #                         a_x_r_dup = return_duplicates(a_x_r)
    #                         if len(c.dup) > len(a_x_r_dup):
    #                             c2.aligned_array = a_x_r
    #                             a_x = a_x_r
    #                             c2.dup = a_x_r_dup
    #                         else:
    #                             if c1.is_terminal():
    #                                 continue
    #                             else:
    #                                 diff = len(c.dup) - len(c1.dup)
    #                                 if diff > 0:
    #                                     a_x_r = self.reversed_alignment(c1)
    #                                     a_x_r, a_y_r, _ = self.nw_alignment(c2.aligned_array, a_x_r)
    #
    #                                     for i, val in enumerate(a_x_r):
    #                                         if val in self.skip_letters:
    #                                             a_x_r[i] = a_y_r[i]
    #                                     a_x_r_dup = return_duplicates(a_x_r)
    #                                     if len(c.dup) > len(a_x_r_dup):
    #                                         c1.aligned_array = a_x_r
    #                                         a_x = a_x_r
    #                                         c1.dup = a_x_r_dup
    #
    #             c.aligned_array = a_x
    #
    #         a_x = tree.root.aligned_array
    #     return a_x

    def reversed_alignment(self, clade):
        clades = copy.deepcopy(clade.clades)
        c1 = clades.pop(-1)

        a_x = c1.aligned_array
        while len(clades) > 0:
            c2 = clades.pop(-1)
            a_x, a_y, _ = self.nw_alignment(c2.aligned_array, a_x)

            for i, val in enumerate(a_x):
                if val in self.skip_letters:
                    a_x[i] = a_y[i]
        return a_x

    def grid_alignment(self, tree=None):
        full_align = self.full_alignment(tree=tree)
        aligned_arrays = pd.DataFrame(columns=list(reversed(range(1, len(full_align) + 1))))
        # print(len(full_align))
        # print(full_align)
        for i, x in enumerate(self.arrays):
            _, a_x, _ = self.nw_alignment(full_align, x)
            # print(a_x)
            aligned_arrays.loc[self.names[i]] = a_x
        return aligned_arrays, full_align


class GAPNWMSA(NWMSA):
    def __init__(self, arrays, names, match=1000, gap=-.1, mismatch=-10000, gap_cont=0, skip_letters=None):
        """
        NOTE: 'Breakpoint' can NOT return actual alignments at the moment (will return an error), just the distance!
        Actually affine gap nw algorithm, probably similar to breakpoint but surely not the same.
        :param arrays:
        :param names:
        :param match:
        :param gap:
        :param mismatch:
        :param gap_cont:
        :param skip_letters:
        """
        super().__init__(arrays, names, match=match, gap=gap, mismatch=mismatch, skip_letters=skip_letters)
        self.gap_cont = gap_cont

    def gap_fct(self, k):
        if k == 0:
            f = self.gap + self.gap_cont
        else:
            f = self.gap_cont
        return f

    def nw_score(self, x, y):
        a, b, c = self.nw_matrix(x, y)
        return max(a[-1, -1], b[-1, -1], c[-1, -1])

    def nw_matrix(self, x, y):
        a = np.zeros((len(x) + 1, len(y) + 1))
        b = np.zeros((len(x) + 1, len(y) + 1))
        c = np.zeros((len(x) + 1, len(y) + 1))
        for i in range(1, len(x) + 1):
            b[i, 0] = -np.infty
            c[i, 0] = -np.infty
            a[i, 0] = a[i - 1, 0] + self.gap_fct(i - 1)
        for j in range(1, len(y) + 1):
            c[0, j] = -np.infty
            b[0, j] = -np.infty
            a[0, j] = a[0, j - 1] + self.gap_fct(j - 1)

        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                b[i, j] = max(a[i - 1, j] + self.gap_fct(0),
                              b[i - 1, j] + self.gap_fct(1))
                c[i, j] = max(a[i, j - 1] + self.gap_fct(0),
                              c[i, j - 1] + self.gap_fct(1))
                a[i, j] = max(a[i - 1, j - 1] + self.match_mismatch(x[i - 1], y[j - 1]),
                              b[i, j],
                              c[i, j])
        return a, b, c

    def match_mismatch(self, a, b):
        if a in self.skip_letters and b in self.skip_letters:
            return 0
        if a == b:
            return self.match
        else:
            return self.mismatch


class JaccardSimilarity(NWMSA):
    """
    Uses actual pairwise breakpoint distance, i.e. number of differing neighbors of equal spacers.
    """

    def __init__(self, arrays, names, match=10, gap=-0.1, mismatch=-10000, skip_letters=None):
        super().__init__(arrays, names, match=match, gap=gap, mismatch=mismatch, skip_letters=skip_letters)

    def _pairwise(self, seq1, seq2):
        # normieren?
        # print(self.nw_score(seq1, seq2))
        return 1 - self.nw_score(seq1, seq2)

    def nw_score(self, x, y):
        inter = 0
        union = 0
        for x_i, y_i in zip(x, y):
            if x_i in self.skip_letters and y_i in self.skip_letters:
                continue
            if x_i == y_i:
                inter += 1
                union += 1
            else:
                union += 1
        if union == 0:
            return 1
        return inter / union


class RealRealBreakpoint(NWMSA):
    """
    Uses actual pairwise breakpoint distance, i.e. number of differing neighbors of equal spacers.
    """

    def __init__(self, match=10, gap=-0.1, mismatch=-10000, skip_letters=None):
        super().__init__(match=match, gap=gap, mismatch=mismatch, skip_letters=skip_letters)

    def _pairwise(self, seq1, seq2):
        # normieren?
        return self.nw_score(seq1, seq2)

    def nw_score(self, x, y):
        score = 0
        clean_x = [s for s in x if s not in self.skip_letters]
        clean_y = [s for s in y if s not in self.skip_letters]
        if len(clean_x) == 0 and len(clean_y) == 0:
            return 0
        equal_spacer_pos, following_spacer, prev_spacer_feq_spacer = self.get_equal_spacer_pos(clean_x, clean_y)
        if equal_spacer_pos:
            for neighbors in following_spacer:
                n_x, n_y = neighbors
                if n_x != n_y:
                    score += 1
            n_x, n_y = prev_spacer_feq_spacer
            if n_x != n_y:
                score += 1
        else:
            return 2
        return score

    def get_equal_spacer_pos(self, x, y):
        equal_spacer_pos = []
        following_spacer = []
        prev_spacer_feq_spacer = None
        for i, x_i in enumerate(x):
            for j, y_j in enumerate(y):
                if x_i == y_j:
                    equal_spacer_pos.append([i, j])
                    if i + 1 >= len(x):
                        fs_x = None
                    else:
                        fs_x = x[i + 1]
                    if j + 1 >= len(y):
                        fs_y = None
                    else:
                        fs_y = y[j + 1]
                    following_spacer.append([fs_x, fs_y])
        if equal_spacer_pos:
            if equal_spacer_pos[0][0] <= 0:
                fs_x = None
            else:
                fs_x = x[equal_spacer_pos[0][0] - 1]
            if equal_spacer_pos[0][1] <= 0:
                fs_y = None
            else:
                fs_y = y[equal_spacer_pos[0][1] - 1]
            prev_spacer_feq_spacer = [fs_x, fs_y]

        return equal_spacer_pos, following_spacer, prev_spacer_feq_spacer


class Breakpoint(NWMSA):
    """
    Works with alignment counts the mutually exclusive gaps. Should be similar to RealRealBreakpoint for
    """

    def __init__(self, arrays, names, match=10, gap=-0.1, mismatch=-10000, skip_letters=None):
        super().__init__(arrays, names, match=match, gap=gap, mismatch=mismatch, skip_letters=skip_letters)

    def _pairwise(self, seq1, seq2):
        # normieren?
        return self.nw_score(seq1, seq2)

    def nw_score(self, x, y):
        if len(x) != len(y):
            raise Exception(
                'len(x)={} is not equal to len(y)={}, i.e. there is no full alignment.'.format(len(x), len(y)))
        score = 0
        # 0 gap in x, 1 gap in y, 2 gap both no gap, 3 gap both
        if len(x) > 0:
            gap = self.gap_fct(x[0], y[0])
            if gap < 2:
                score += 1
        else:
            return 0
        ls_gap = []
        for x_i, y_i in zip(x[1:], y[1:]):
            gap_next = self.gap_fct(x_i, y_i)
            if gap_next == 2 and x_i != y_i:
                score += 1
            # print('gap', gap)
            # print('gapnext', gap_next)
            if gap != gap_next and gap_next < 2:
                score += 1
            # print(gap)
            # print('score', score)
            gap = gap_next
            ls_gap.append(gap)
        # set_gap = set(ls_gap)
        # if len(set_gap) == 1:
        #     if set_gap.pop() < 2:
        #         score = 1
        #     else:
        #         score = 0
        return score

    def gap_fct(self, x_i, y_i):
        if x_i in self.skip_letters:
            if y_i in self.skip_letters:
                gap = 3
            else:
                gap = 0
        elif y_i in self.skip_letters:
            gap = 1
        else:
            gap = 2
        return gap
        # clean_x = [s for s in x if s not in self.skip_letters]
        # clean_y = [s for s in y if s not in self.skip_letters]
        # score = 0
        # first = True
        # k = 0
        # for i in range(len(clean_x)):
        #     val = clean_x[i]
        #     j = k
        #     while j < len(clean_y):
        #         if val == clean_y[j]:
        #             if first:
        #                 if (j == 0) or (i == 0):
        #                     first = False
        #                     if not (j == 0 and i == 0):
        #                         score += 1
        #                 elif clean_y[j - 1] != clean_x[i - 1]:
        #                     first = False
        #                     score += 1
        #             if (j + 1 >= len(clean_y)) or (i + 1 >= len(clean_x)):
        #                 if not (j + 1 >= len(clean_y) and i + 1 >= len(clean_x)):
        #                     score += 1
        #             elif clean_y[j + 1] != clean_x[i + 1]:
        #                 score += 1
        #             k = j + 1
        #         j += 1
        #
        # if first:
        #     score = len(clean_x) + len(clean_y)
        # return score
        # print('clean x', clean_x)
        # print('clean y', clean_y)
        # adj_x = []
        # adj_y = []
        # for i in range(len(clean_x) - 1):
        #     adj_x.append([clean_x[i], clean_x[i + 1]])
        # for j in range(len(clean_y) - 1):
        #     pair = [clean_y[j], clean_y[j + 1]]
        #     idx_to_del = []
        #     for i in range(len(adj_x)):
        #         if pair == adj_x[i]:
        #             idx_to_del.append(i)
        #     if idx_to_del:
        #         for idx in idx_to_del:
        #             del adj_x[idx]
        #     else:
        #         adj_y.append(pair)
        # print('res2', adj_x, adj_y)
        # print('res', len(adj_x + adj_y))
        # return len(adj_x) + len(adj_y)


def jaccard_distance_2(a, b, order):
    return order


def jaccard_distance(a, b):
    cap = len(set(a) - set(b))
    cup = len(set(a) + set(b))
    return cap / cup


def breakpoint_distance(a, b):
    dist = 0
    return dist


def get_distance_matrix(arrays, distance=None):
    return 0


def return_duplicates(array):
    series = pd.Series(array)
    mask = series.duplicated(keep='last')
    dups = series[mask]
    return dups
