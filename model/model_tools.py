import numpy as np
import copy
import scipy.optimize
import scipy.stats
import pandas as pd
import collections
import os

from model.helpers import misc
from model.spacer_visualization import visualization as vis

EPS = 1e-40
BIG_EPS = 1e20


def minimize_fct(lh_fct, start_values, bounds=None, method='L-BFGS-B'):
    opt = scipy.optimize.minimize(lh_fct, start_values, bounds=bounds, method=method)
    return opt


def minimize_scalar_fct(lh_fct, bounds=None, method='bounded'):
    return scipy.optimize.minimize_scalar(lh_fct, bounds=bounds, method=method)


def compute_lh_ratio_of_multiple_trees(ls_data, method='Nelder-Mead', filter_by=False, lambdifyed_lh_fct=None):
    lh_fct_0 = lambda x: -sum([compute_tree_lh(x, 1.0, d[0], d[1], d[2], log=True) for d in ls_data])

    if lambdifyed_lh_fct is None:
        lh_fct_1 = lambda x: -sum([compute_tree_lh(x[0], x[1], d[0], d[1], d[2], log=True)
                                   for d in ls_data])
    else:
        lh_fct_1 = lambda x: -sum([compute_tree_lh_for_given_lh_fct(x[0], x[1], d[0], d[1], d[2], log=True,
                                                                    lambdifyed_lh_fct=lambdifyed_lh_fct)
                                   for d in ls_data])

    start_values_1 = np.array([0.5, 1.5])
    bnds_0 = (0.0, 10000)
    bnds_1 = ((0.0, None), (1.0, None))
    # print(len(start_values_0), len(bnds_0))

    max_0 = minimize_scalar_fct(lh_fct_0, bounds=bnds_0, method='bounded')
    max_1 = minimize_fct(lh_fct_1, start_values_1, bounds=bnds_1, method=method)

    # lh_0 = compute_tree_lh(max_0.x, 1, ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses, log=False)
    # lh_1 = compute_tree_lh(max_1.x[0], max_1.x[1], ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses, log=False)
    lh_0 = np.exp(-max_0.fun)
    lh_1 = np.exp(-max_1.fun)

    return lh_0, lh_1, max_0, max_1


def compute_tree_lh(loss_rate, alpha, ls_bl, nb_survivors, nb_max_length_losses, log=True,
                    rec_tree=None, gains=None):
    tree_lh = []
    if rec_tree is None:
        for i, (bl, n, mll) in enumerate(zip(ls_bl, nb_survivors, nb_max_length_losses)):
            loss_prob = max(1 - np.exp(- loss_rate * bl), 0)
            max_mll = max(mll) if mll else 0
            ls_f = eff_mult_loss_lh(loss_prob, alpha, max_mll, log=True)
            ls_f[np.isnan(ls_f) | np.isinf(ls_f)] = -BIG_EPS
            # print('lossrate', loss_rate)
            # print('bl', bl)
            # print('lp', loss_prob)
            # print('bl', bl, n, mll)
            ln_keep = - loss_rate * bl * n
            # print('ln_keep', ln_keep)
            # print('lsf', ls_f)
            # print('mll', mll)
            # print('f', ls_f)
            ls_ln_mult_loss_lh = [ls_f[k_i - 1] for k_i in mll]

            ln_mult_loss_lh = sum(ls_ln_mult_loss_lh)
            # print('fki', ls_ln_mult_loss_lh)
            # print('final lh before adding keep', ln_mult_loss_lh)
            tree_lh.append(ln_mult_loss_lh + ln_keep)
            # print('final', ln_mult_loss_lh + ln_keep)
    else:
        # print('gains', gains)
        # print('gains', gains[1:])
        eff_lr = loss_rate
        ls_tree_ln_u_lh = rec_tree.compute_unobserved_lh(eff_lr)
        for i, (bl, n, mll, g) in enumerate(zip(ls_bl, nb_survivors, nb_max_length_losses, gains[1:])):
            loss_prob = max(1 - np.exp(- loss_rate * bl), 0)
            max_mll = max(mll) if mll else 0
            ls_f = eff_mult_loss_lh(loss_prob, alpha, max_mll, log=True)

            ls_f[np.isnan(ls_f) | np.isinf(ls_f)] = -BIG_EPS
            ln_keep = - loss_rate * bl * n
            ls_ln_mult_loss_lh = [ls_f[k_i - 1] for k_i in mll]

            ln_u_lh = ls_tree_ln_u_lh[i + 1]
            if loss_rate != 0 and ln_u_lh is not None:
                ln_u_lh = log1mexp(ls_tree_ln_u_lh[i + 1])
                ls_ln_unobserved_lh_inv = [ln_u_lh] * g
                ls_ln_mult_loss_lh += [-lh for lh in ls_ln_unobserved_lh_inv]
            ln_mult_loss_lh = sum(ls_ln_mult_loss_lh)
            tree_lh.append(ln_mult_loss_lh + ln_keep)
        # print('tree_lh', sum(tree_lh))
        root_gain = gains[0]
        if loss_rate != 0:
            ln_u_lh = log1mexp(ls_tree_ln_u_lh[0])
            ls_ln_unobserved_lh_inv = [ln_u_lh] * root_gain
            tree_lh.append(sum([-lh for lh in ls_ln_unobserved_lh_inv]))
    return sum(tree_lh) if log else np.exp(sum(tree_lh))


def log1mexp(x):
    """
    Computes log(1 - exp(x)) in a numerically stable way. Note, this is only viable for abs(x) < 1 and x < 0.
    log(2) \approx 0.693 seems to be the optimal cutoff (for highest precision).
    :param x:
    :return:
    """
    # print('log1mexp1', np.log1p(-np.exp(x)))
    # print('log1mexp2', np.log(-np.expm1(x)))
    if x < -0.693:
        return np.log1p(-np.exp(x))
    else:
        return np.log(-np.expm1(x))


def compute_tree_lh_for_given_lh_fct(loss_rate, alpha, ls_bl, nb_survivors, nb_max_length_losses, log=True,
                                     lambdifyed_lh_fct=None):
    """
    Computes blm-lh of a tree.
    :param lambdifyed_lh_fct:
    :param loss_rate:
    :param alpha:
    :param ls_bl:
    :param nb_survivors:
    :param nb_max_length_losses:
    :param log:
    :return:
    """
    tree_lh = []
    for bl, n, mll in zip(ls_bl, nb_survivors, nb_max_length_losses):
        # loss_prob = max(1 - np.exp(- loss_rate * bl), 0)
        max_mll = max(mll) if mll else 0

        # Could remove RuntimeWarning by
        # res = np.log(m, out=np.zeros_like(m), where=(m!=0))
        # or
        # with errstate(divide='ignore'):
        #     res = np.log(m)
        ls_f = np.log([max(f(bl, loss_rate, alpha), EPS) for f in lambdifyed_lh_fct[:max_mll + 1]])
        ls_f[np.isnan(ls_f) | np.isinf(ls_f)] = -BIG_EPS
        ln_keep = - loss_rate * bl * n

        # Handle blocks that are longer than the available likelihood functions. Divide them in max length blocks
        # + Rest.
        if max_mll + 1 > len(ls_f):
            ls_ln_mult_loss_lh = []
            for k_i in mll:
                val = 0
                val += ls_f[k_i % (len(ls_f) - 1)]
                incl_number = k_i // (len(ls_f) - 1)
                # print('incl_number', incl_number, 'k_i % len(ls_f)', k_i % (len(ls_f) - 1))
                for j in range(incl_number):
                    val = incl_number * ls_f[-1]
                ls_ln_mult_loss_lh.append(val)
        else:
            ls_ln_mult_loss_lh = [ls_f[k_i] for k_i in mll]
        ln_mult_loss_lh = sum(ls_ln_mult_loss_lh)
        # print('fki', ls_ln_mult_loss_lh)
        # print('final lh before adding keep', ln_mult_loss_lh)
        tree_lh.append(ln_mult_loss_lh + ln_keep)
        # print('final', ln_mult_loss_lh + ln_keep)
    return sum(tree_lh) if log else np.exp(sum(tree_lh))


def eff_mult_loss_lh(loss_prob, alpha, k_i, log=True):
    geom_k = [geom(alpha, k) for k in range(1, k_i + 1)]
    # print(alpha)
    lh = mult_loss_lh(loss_prob, alpha, geom_k)
    lh = [max(EPS, x) for x in lh]
    # if lh > 1 or lh <= 0:
    #     print('alpha', alpha)
    #     print('prob_k', prob_k)
    #     print('lh', lh)
    return np.log(lh) if log else lh


def mult_loss_lh(loss_prob, alpha, geom_k, log=True):
    """
    Reproduced formula:
    f(k, \rho, \alpha) = P(loss)*(sum_{j=1}^{k} (k + 1 - j)*f(k-j)*g_\alpha (j)
                       = P(loss)*((sum_{j=1}^{k-1} (k + 1 - j)*f(k-j)*g_\alpha (j)
                       + 1 * f(0)(=1) * g_\alpha (k))
    The Term (k + 1 - j) seems to be wrong, since lh increases above 1 in few steps. f diverges for k \to \infty, i.e.
    (k + 1 - j) *... increase faster than the probabilities (P(loss), g_\alpha (k)).
    See \alpha = 1: Then f(k) = P(loss)^k * \prod_{j=1}^k j (= k!), i.e. f(k) < 1 only if k, P(loss) small enough.
    In general need condition P(loss)^k * g_\alpha(something) * ... * g_\alpha(something) is smaller than some product
    of (k + 1 - j)..., s.t. f(k) < 1.
    :param loss_prob:
    :param alpha:
    :param geom_k:
    :param log:
    :return:
    """
    k = len(geom_k)
    if k == 0:
        return []
    # if k == 1:
    #     return [loss_prob * geom_k[-1]]
    # print('alpha', alpha)
    # print('loss_prob', loss_prob)
    # print('prob_K', geom_k)
    # for j in range(0, k - 1):
    #     val = geom_k[j] * mult_loss_lh(loss_prob, alpha, geom_k[:-(j+1)])
    #     ls.append(val)
    # ls = [prob_k[j-1] * mult_loss_lh(loss_prob, alpha, prob_k[:-j]) for j in range(1, k)]
    f = mult_loss_lh(loss_prob, alpha, geom_k[:-1])
    ls = [geom_k[j] * f[-(j + 1)] for j in range(0, k - 1)]
    # new f computation with (k + 1 - j)/(k - (j - 1) factor, have to be careful with indices and stuff
    # ls = [geom_k[j] * f[-(j + 1)] * (k - j) for j in range(0, k - 1)]
    f.append(loss_prob * (geom_k[-1] + sum(ls)))
    # print('prob_k-1', geom_k[-1])
    # print('ls', ls)
    # print('sum', sum(ls))
    # print('ret', ret)
    return f


def geom(alpha, k):
    if alpha == 1:
        if k == 1:
            return 1
        else:
            return 0
    p = 1 / alpha
    return (1 - p) ** (k - 1) * p


def create_df_alignment(ls_arrays, ls_names, topological_order=None):
    """
    I hate this function, should remove it, because it isn't beneficial except for printing the alignment.
    But I'm also using dfs in other functions so...
    :param ls_arrays:
    :param ls_names:
    :param topological_order:
    :return:
    """
    if topological_order:
        aligned_arrays = pd.DataFrame(columns=topological_order)
        for i, x in enumerate(ls_arrays):
            aligned_array = [s if s in x else '-' for s in topological_order]
            aligned_arrays.loc[ls_names[i]] = aligned_array
    else:
        # Creates single 'empty' column as placeholder, hopefully fixes issues with empty dataframes...
        if len(ls_arrays[0]) == 0:
            raise ValueError('The spacer alignment is empty!')
        else:
            aligned_arrays = pd.DataFrame(columns=list(range(1, len(ls_arrays[0]) + 1)))
            for i, x in enumerate(ls_arrays):
                # print(len(x))
                # print(len(aligned_arrays.columns))
                aligned_arrays.loc[ls_names[i]] = x
    return aligned_arrays


def remove_unobserved_spacer(df_aligned, tree, gain_dict, loss_dict, gain_loss_event_dict=None, logger=None):
    """
    Could also only search for only '-' columns in alignment, should be cheaper.
    :param tree:
    :param gain_loss_event_dict:
    :param df_aligned:
    :param top_order:
    :param gain_dict:
    :param loss_dict:
    :return:
    """
    top_order = []
    top_order_col = []
    for col in df_aligned.columns:
        unique_val = next((item for item in df_aligned[col] if item != '-'), None)
        if unique_val is not None:
            top_order.append(unique_val)
            top_order_col.append(col)

    # Determines the nodes that are considered to observe spacers, i.e. the leafs in most cases.
    considered_nodes = [node.name for node in tree.get_terminals()]
    observed_spacers = set()
    for name in considered_nodes:
        observed_spacers = observed_spacers.union(set(df_aligned.loc[name]))
    observed_spacers = observed_spacers - {'-'}

    new_top_order = []
    col_to_drop = []
    for s, col in zip(top_order, top_order_col):
        if s not in observed_spacers:
            col_to_drop.append(col)
        else:
            new_top_order.append(s)
    new_df = df_aligned.drop(col_to_drop, axis=1)
    if logger:
        logger.info(f'Original number of spacers: {df_aligned.shape[1]} ; '
                    f'Number of spacers w/o unobserved: {new_df.shape[1]}')
    else:
        print('Original number of spacers: ', df_aligned.shape[1],
              'Number of spacers w/o unobserved: ', new_df.shape[1])
    new_gain_dict = {}
    new_loss_dict = {}
    for key in gain_dict.keys():
        new_gain_dict[key] = [gain for gain in gain_dict[key] if gain in observed_spacers]
        new_loss_dict[key] = [loss for loss in loss_dict[key] if loss in observed_spacers]

    if gain_loss_event_dict:
        for key, ls_events in gain_loss_event_dict.items():
            for i, event in enumerate(ls_events):
                event[2] = [s for s in event[2] if s in observed_spacers]
                if not event[2]:
                    del ls_events[i]
                event[4] = [s for s in event[4] if s in observed_spacers]

    return new_df, new_top_order, new_gain_dict, new_loss_dict


def find_wrong_order_gains_dict(rec_gains, tree, save_to_tree=True, given_upgraph_order=None):
    """
    Compares previously gained spacers to known younger spacers (upgraph order) of each gain at current node. If there
    is an intersection --> contradiction.
    :param rec_gains:
    :param tree:
    :param save_to_tree:
    :param given_upgraph_order:
    :return:
    """
    contradictions = {}
    gained_spacers = {}
    for c in tree.find_clades(order='preorder'):
        if c.up is None:
            gained_spacers[c.name] = set(rec_gains[c.name])
            continue

        for gain in rec_gains[c.name]:
            if len(gained_spacers[c.up.name].intersection(given_upgraph_order[gain])) > 0:
                contradictions[c.name] = [gain] + contradictions.get(c.name, [])
        gained_spacers[c.name] = gained_spacers[c.up.name].union(set(rec_gains[c.name]))
        if save_to_tree:
            c.contra_gain = contradictions.get(c.name, [])
    return contradictions


# def detect_loop(leaf_arrays):
#     unique_spacers_pos = {spacer: pos for pos, spacer in enumerate(list({x for ls in leaf_arrays for x in ls}))}
#     adj_matrix = np.zeros((len(unique_spacers_pos), len(unique_spacers_pos)))
#     # generate adjacency matrix
#     for ls in leaf_arrays:
#         for i, spacer in enumerate(ls):
#             if i == 0:
#                 continue
#             adj_matrix[unique_spacers_pos[ls[i - 1]], unique_spacers_pos[ls[i]]] = 1
#     # plot graph
#     label_dict = {pos: spacer for spacer, pos in unique_spacers_pos.items()}
#
#     g = nx.from_numpy_array(adj_matrix, create_using=nx.DiGraph)
#     cycles = nx.simple_cycles(g)
#     labeled_cycles = []
#     for c in cycles:
#         s_c = []
#         for n in c:
#             s_c.append(label_dict[n])
#         labeled_cycles.append(s_c)
#     return labeled_cycles


def determine_order_dict(leaf_arrays, save_path=None, plot_order=True, exclude_self=True):
    """
    Creates a partial spacer insertion order graph from leaf arrays and plots it. Then determines the order of
    spacers: determines the subgraph, older spacers, and upgraph, younger spacers, for each spacer and
    saves them to respective dict[spacer]. Also determines a topological order. Initiates visualization of order.
    :param leaf_arrays:
    :param save_path:
    :param plot_order:
    :param exclude_self: If True, excludes the spacer from sets in order dicts. I.e. dict[spacer] = set() instead of
    dict[spacer] = {spacer} for oldest/youngest spacer in subgraph/upgraph dict, respectively.
    :return:
    """
    unique_spacers_pos = {spacer: pos for pos, spacer in enumerate(list({x for ls in leaf_arrays for x in ls}))}
    adj_matrix = np.zeros((len(unique_spacers_pos), len(unique_spacers_pos)))
    # generate adjacency matrix
    for ls in leaf_arrays:
        for i, spacer in enumerate(ls):
            if i == 0:
                continue
            adj_matrix[unique_spacers_pos[ls[i - 1]], unique_spacers_pos[ls[i]]] = 1
    # plot graph
    label_dict = {pos: spacer for spacer, pos in unique_spacers_pos.items()}

    # find spacers/nodes with no incoming edges (roots)
    s = []
    for j in range(adj_matrix.shape[1]):
        if all(a == 0 for a in adj_matrix[:, j]):
            s += [j]

    # get one topological order
    adj_matrix_ = copy.deepcopy(adj_matrix)
    s_ = copy.deepcopy(s)
    topological_order = misc.topological_sort(s_, adj_matrix_, label_dict)

    # right sets (subtrees):
    subgraph_order_dict = {pos: set() for pos in label_dict.keys()}
    for root in s:
        subgraph_order_dict[root] = get_subgraph_sets(root, adj_matrix, subgraph_order_dict)

    if exclude_self:
        for key, subgraph in subgraph_order_dict.items():
            subgraph_order_dict[key] = subgraph - {key}
    # return names
    named_subgraph_order_dict = {}
    level_dict = {}
    for key, subgraph in subgraph_order_dict.items():
        named_subgraph_order_dict[label_dict[key]] = set(map(lambda x: label_dict[x], subgraph))
        level_dict[key] = -len(subgraph)

    # left sets:
    upgraph_order_dict = {}
    for spacer in topological_order:
        spacer_left_set = set()
        for key, subgraph in named_subgraph_order_dict.items():
            if spacer in subgraph:
                spacer_left_set.add(key)
        upgraph_order_dict[spacer] = spacer_left_set

    if plot_order:
        path = os.path.split(save_path)[0]
        path_2 = os.path.split(save_path)[1]
        vis.plot_order_w_graphviz(adj_matrix, label_dict=label_dict, do_show=False,
                                  save_folder=path,
                                  graph_name=path_2.split('.')[0],
                                  file_name=path_2 + '.dot',
                                  color_dict=None)

    return (named_subgraph_order_dict, upgraph_order_dict), topological_order, (adj_matrix, label_dict)


def get_subgraph_sets(n, adj_m, save_dict):
    """
    Do I want to have <= or < spacers in subtree? Currently, <= for all. Exclusion is manually facilitated in
    determine_order_dict (exclude_self).
    :param n:
    :param adj_m:
    :param save_dict:
    :return:
    """
    if save_dict[n]:
        return save_dict[n]
    elif all(j == 0 for j in adj_m[n, :]):
        save_dict[n] = {n}
        return {n}
    else:
        total_subgraph = {n}
        for col, j in enumerate(adj_m[n, :]):
            if j == 1:
                subgraph = get_subgraph_sets(col, adj_m, save_dict)

                total_subgraph = total_subgraph.union(subgraph)
                # print(total_subtrees)
        save_dict[n] = total_subgraph
        return total_subgraph


def get_fix_mismatches_df_alignment(df, reverse_walkthrough=False):
    """
    Fixes errors produced by mafft. There are columns with multiple spacers, that need to be separated.
    More problematic are artificial duplications produced by mafft, i.e. columns with the same name that actually can
    be combined, if the PSIO is considered.
    :param df:
    :return:
    """
    if df.empty:
        return ValueError('The spacer alignment is empty!')
    mismatches = []
    count_mismatches = 0

    ls_resolved_mismatches_col = []
    series_name = 0
    old_uniques = set()
    set_resolved_mismatch_names = set()
    # Separate columns with multiple spacers. They are moved accordingly to what is to the left or right of them, to not
    # introduce additional artificial duplicates.

    # if reverse_walkthrough:
    #     df = df[df.columns[::-1]]

    for idx in range(df.shape[1]):
        series = df.iloc[:, idx]
        uniques = set(series) - {'-'}

        next_uniques = set(df.iloc[:, idx + 1]) - {'-'} if idx + 1 < df.shape[1] else set()
        if len(uniques) > 1:
            inter_w_old = uniques.intersection(old_uniques)
            inter_w_next = uniques.intersection(next_uniques)
            inter_all = inter_w_old.intersection(inter_w_next)
            remaining_uniques = uniques - inter_w_old - inter_w_next
            inter_w_old -= inter_all
            inter_w_next -= inter_all
            ls_first = []
            ls_last = []
            ls_middle = []

            for s in remaining_uniques:
                resolved_mismatches_col = pd.Series([s if series.iloc[i] == s else '-' for i in range(len(series))],
                                                    name=series_name, index=series.index)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_last.append(resolved_mismatches_col)
            for s in inter_all:
                resolved_mismatches_col = pd.Series([s if series[i] == s else '-' for i in range(len(series))],
                                                    name=series_name, index=series.index)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_middle.append(resolved_mismatches_col)
            for s in inter_w_next:
                resolved_mismatches_col = pd.Series([s if series[i] == s else '-' for i in range(len(series))],
                                                    name=series_name, index=series.index)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_last.append(resolved_mismatches_col)
            for s in inter_w_old:
                resolved_mismatches_col = pd.Series([s if series[i] == s else '-' for i in range(len(series))],
                                                    name=series_name, index=series.index)
                set_resolved_mismatch_names.add(s)
                series_name += 1
                ls_first.append(resolved_mismatches_col)
            ls_resolved_mismatches_col += ls_first + ls_middle + ls_last

        else:
            series.name = series_name
            series_name += 1
            ls_resolved_mismatches_col.append(series)

        old_uniques = set(ls_resolved_mismatches_col[-1]) - {'-'}
    resolved_mismatches_df = pd.concat(ls_resolved_mismatches_col, axis=1)

    # if reverse_walkthrough:
    #     resolved_mismatches_df = resolved_mismatches_df[resolved_mismatches_df.columns[::-1]]

    candidates = [(col, next(iter(set(resolved_mismatches_df[col]) - {'-'})))
                  for col in resolved_mismatches_df.columns if len(set(resolved_mismatches_df[col]) - {'-'}) != 0
                  ]

    duplicates = [item for item, count in collections.Counter([val[1]
                                                               for val in candidates]).items()
                  if count > 1]

    # Check if duplicates can be joined by looking at the order between them.
    dict_of_ls_cols_to_join = dict()
    min_idx = resolved_mismatches_df.shape[1]
    max_idx = 0
    # min_col, max_col = None, None
    for dup in duplicates:
        columns_not_to_join = []
        all_cols = [(i, col) for i, (col, val) in enumerate(candidates) if val == dup]
        new_min_idx, new_min_col = min(all_cols, key=lambda x: x[0])
        new_max_idx, new_max_col = max(all_cols, key=lambda x: x[0])
        if new_min_idx < min_idx:
            min_idx = new_min_idx
            # min_col = list(resolved_mismatches_df.columns)[min_idx]
        if new_max_idx > max_idx:
            max_idx = new_max_idx
            # max_col = list(resolved_mismatches_df.columns)[max_idx]
        in_between_cols = copy.deepcopy(resolved_mismatches_df.iloc[:, min_idx:max_idx + 1])
        ls_renamed_dup = dict()
        columns = set(in_between_cols.columns)
        for d in duplicates:
            # print('d', d)
            dup_counter = 0
            # columns that can be joined w/o issue
            joinables = dict_of_ls_cols_to_join.get(d, [])
            # columns that can be joined w/o issue and are actually in between columns,
            same_name = []
            # new names for same_name spacers
            ls_dup_name = []
            i = 0
            for ls in joinables:
                same_name.append(set([a for a in ls if a in columns]))
                ls_dup_name += [i]
                i += 1
            # Give joinable spacers the same name
            if same_name:
                for j, s in enumerate(same_name):
                    for col in s:
                        for idx in in_between_cols.index:
                            val = in_between_cols.loc[idx, col]
                            in_between_cols.loc[idx, col] = d + '_' + str(ls_dup_name[j]) if val == d else val
            else:
                # Give not joinables individual names
                for col in in_between_cols.columns:
                    # if col not in skip_renaming and d in set(in_between_cols[col]):
                    if d in set(in_between_cols[col]):
                        # print(col)
                        for idx in in_between_cols.index:
                            val = in_between_cols.loc[idx, col]
                            in_between_cols.loc[idx, col] = d + '_' + str(dup_counter) if val == d else val
                        # Remember duplications that need be transformed back.
                        if d == dup:
                            ls_renamed_dup[d + '_' + str(dup_counter)] = col
                        dup_counter += 1
        ls_aligned_arrays = [[value for value in ls if value != '-'] for ls in in_between_cols.values]
        # try:
        (named_subgraph_order_dict, upgraph_order_dict), _, _ = determine_order_dict(ls_aligned_arrays,
                                                                                     save_path=os.path.join(
                                                                                         '../0_result_folder',
                                                                                         'testing'),
                                                                                     plot_order=False,
                                                                                     )

        # Find if up and down order of one individual overlap with any other duplicate candidate,
        # if so they are divided and can't be joined!
        # Should implement this with permutation invariance, i.e. only one of (1, 2) (2, 1).
        for rd, col in ls_renamed_dup.items():
            set_wo_rd = set(ls_renamed_dup.keys()) - {rd}
            subgraph_order_rd = named_subgraph_order_dict[rd]
            upgraph_order_rd = upgraph_order_dict[rd]
            inter_subgraph = set_wo_rd.intersection(subgraph_order_rd)
            inter_upgraph = set_wo_rd.intersection(upgraph_order_rd)
            set_d_not_to_join = inter_upgraph.union(inter_subgraph)
            col_not_to_join = [ls_renamed_dup[s] for s in set_d_not_to_join]

            # col is not allowed to be joined with c.
            for c in col_not_to_join:
                columns_not_to_join.append((col, c))

        ls_usable_cols = list(ls_renamed_dup.values())
        # Put join candidates together for maximum length connected sets, i.e.
        # don't allow not joinable pairs from above
        while ls_usable_cols:
            new_ls = []
            to_keep = []
            for i in range(len(ls_usable_cols)):
                candidate_col = ls_usable_cols[i]
                addable = all([True if (candidate_col, c) not in columns_not_to_join else False for c in new_ls])
                if addable or not new_ls:
                    new_ls.append(candidate_col)
                else:
                    to_keep.append(candidate_col)
            dict_of_ls_cols_to_join[dup] = dict_of_ls_cols_to_join.get(dup, []) + [new_ls]
            ls_usable_cols = to_keep

    top_order = []
    for col in resolved_mismatches_df.columns:
        unique_col = pd.unique([a for a in resolved_mismatches_df[col] if a != '-'])
        if len(unique_col) > 0:
            if len(unique_col) > 1:
                raise Exception('You fucked up! Mismatch in column' + str(col))
            else:
                val = list(unique_col)
        else:
            val = []
        top_order += val

    mask = [True]
    series = pd.Series(top_order)
    dups = []
    count_dups = 0
    while any(mask):
        mask = series.duplicated(keep='last')
        series = series[mask]
        unique_series = pd.unique(series)
        dups += (list(unique_series))
        count_dups += len(unique_series)
    return resolved_mismatches_df, top_order, dict_of_ls_cols_to_join, (count_dups, count_mismatches, mismatches, dups)


def differences_child_parent(child_profile, parent_profile):
    diff_p_ch = child_profile - parent_profile
    # should be careful this relies on maximal length of the profile (for the labeling in the dictionary)
    len_differences = len(diff_p_ch)
    gains = [(len_differences - idx) for idx, val in enumerate(diff_p_ch) if val == 1]
    losses = [(len_differences - idx) for idx, val in enumerate(diff_p_ch) if val == -1]
    gains = list(reversed(gains))
    losses = list(reversed(losses))
    return gains, losses


def differences_child_parent_pos(child_profile, parent_profile):
    diff_p_ch = child_profile - parent_profile
    gains = [idx for idx, val in enumerate(diff_p_ch) if val == 1]
    losses = [idx for idx, val in enumerate(diff_p_ch) if val == -1]
    return gains, losses


def seq2prof(seq, profile_map=None):
    """
    Converts presence absence to likelihood profile (default guaranteed 0 or 1). This goes solely by position.
    1-d -> 2-d array
    :param seq:
    :param profile_map:
    :return:
    """
    if profile_map is None:
        profile_map = {0: [1.0, 0], 1: [0, 1.0]}
    return np.array([profile_map[k] for k in seq])


def prof2seq(profile):
    """
    Converts likelihood profile to presence absence profile. This goes solely by position.
    2-d -> 1-d array
    :param profile:
    :return:
    """
    if len(profile.shape) == 1:
        length = len(profile)
        seq = []
        for i, p in enumerate(profile):
            if p == 1:
                seq.append(length - i)
        seq = np.array(seq)
    else:
        seq = profile.argmax(axis=1).astype(int)
    return seq


def spacer_array_to_binary(spacer_array, alignment=None):
    # if 0 in spacer_array:
    #     raise ValueError('Spacer array contains 0!')
    return np.isin(alignment, test_elements=spacer_array).astype(int)
