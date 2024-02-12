import os
import numpy as np
import pandas as pd
import string
import copy
import collections
import json

from model.data_classes import advanced_tree
from additional_data.additional_scripts.model.evolution import GTR
from model.model_tools import seq2prof, prof2seq, differences_child_parent_pos, \
    compute_tree_lh, \
    minimize_fct, minimize_scalar_fct
from model.spacer_visualization import visualization as vis
from model.helpers import misc
from model import model_tools


class ReconstructionTree(advanced_tree.AdvancedTree):
    def __init__(self, rec_model_config, save_path, root, rooted, model_name=None, lh_fct=None,
                 sim_as_rec=False, sim_gain_loss_dicts=None, hide_unobserved_spacers=False):
        """
        We require the tree to have unique names for each node/branch!
        :param rec_model_config: Only GTR is used and implemented. gain_rates, loss_rates should be chosen such that
        gain_rate is much lower than loss_rate to encourage gaining spacer only a single time.
        :param root:
        :param rooted:
        :param model_name:
        :param lh_fct:
        :param hide_unobserved_spacers: hides unobserved spacers, i.e. spacers that are not found in any leaf. This
        should only occur, if simulations are used.
        """
        super(ReconstructionTree, self).__init__(root=root, rooted=rooted, model_name=model_name)

        self.hide_unobserved_spacers = hide_unobserved_spacers
        self.sim_as_rec = sim_as_rec
        self.sim_gain_loss_dicts = sim_gain_loss_dicts

        self.rec_gain_dict = {}
        self.rec_loss_dict = {}
        self.rec_contra_dict = {}
        self.rec_joint_losses_dict = {}

        self.eps = 1e-40
        self.big_eps = 1e20

        if rec_model_config['model'].lower() == 'gtr':
            self.rec_model = GTR(rec_model_config['gain_rate'], rec_model_config['loss_rate'])
            self.gain_rate = rec_model_config['gain_rate']
            self.loss_rate = rec_model_config['loss_rate']
        else:
            raise NotImplementedError('rec_model_config not supplied or model not implemented!')

        self.spacer_orders = None
        self.spacer_up_orders = None
        self.top_order = None

        self.save_path = save_path

        # Might want to remove this if case
        self.lambdifyed_lh_fct = lh_fct

        self.spacer_names_to_numbers = None
        self.spacer_numbers_to_names = None

        self.duplicated_spacers = None
        # expected to be duplications
        self.rec_duplications_dict = {}
        # expected to be double gains jumps or inversions
        self.rec_rearrangements_dict = {}
        self.rec_reacquisition_dict = {}
        # For now, they are independent acquisitions and a default case.
        self.rec_indep_gain_dict = {}
        self.rec_other_dup_events_dict = {}

        self.rec_dup_rearr_candidates = {}
        self.order_adj_matrix = None

    ####################################################################################################################
    # Get things from the tree.
    def count_moved_spacers(self, prev_gain_dict: dict) -> int:
        """
        Counts the number of spacers where an insertion moved in the tree through the refinement step.
        :param prev_gain_dict:
        :return:
        """
        encountered_spacers = set()
        count_moved_spacers = 0
        for node in self.root.find_clades(order='level'):
            prev_gains = set(prev_gain_dict[node.name])
            curr_gains = set(self.rec_gain_dict[node.name])

            curr_gains_wo_encountered_gains = curr_gains - encountered_spacers
            count_moved_spacers += len(curr_gains_wo_encountered_gains - prev_gains)
            encountered_spacers = encountered_spacers.union(curr_gains)
        return count_moved_spacers

    def get_dict_duplicated_spacers(self):
        """
        Gets the dictionary of duplicated spacers but with spacer names instead of numbers.
        :return:
        """
        return {key: set(map(lambda x: self.spacer_numbers_to_names[x], value))
                for key, value in self.duplicated_spacers.items()}

    def avg_length_leaf_spacers(self):
        """
        Returns the mean length of arrays in the leaves. Only works after reconstruction.
        :return:
        """
        ls_spacers_len = []
        for c in self.root.get_terminals():
            ls_spacers_len.append(sum(c.bin_spacers))
        return np.mean(ls_spacers_len)

    def get_leaf_array_lengths(self):
        """
        Returns a dictionary of all leaf array lengths
        :return:
        """
        dict_lengths = {}
        for leaf in self.root.get_terminals():
            dict_lengths[leaf.name] = len(leaf.spacers)
        return dict_lengths

    def get_data_for_lh_ratio(self, filter_by=False):
        """
        Returns lists of branch lengths, surviving spacers, and lengths of adjacent losses.
        :param filter_by: Allows filtering by the presence of spacers in the root.
        """
        self.generate_possible_losses_sets()
        # Reruns the computation of the maximal length losses to have the right filter.
        self.get_max_length_loss_sets(filter_by=filter_by)
        ls_bl = []
        ls_nb_ex_spacers = []
        ls_nb_max_length_losses = []
        for node in self.root.find_clades(order='preorder'):
            if node.up is None:
                continue
            ls_bl.append(node.branch_length)
            parent_profile = node.up.rec_spacers
            child_profile = node.rec_spacers
            if filter_by:
                # Easier for "contains old spacers" there it is just all present in root.
                root_profile = self.root.rec_spacers
                ls_nb_ex_spacers.append(sum(np.array([1 if (p == 1 and c == 1 and f == 1) else 0 for c, p, f in
                                                      zip(child_profile, parent_profile, root_profile)])))
            else:
                ls_nb_ex_spacers.append(
                    sum(np.array([1 if (p == 1 and c == 1) else 0 for c, p in zip(child_profile, parent_profile)])))
            ls_nb_max_length_losses.append(node.max_length_losses[1])
            # print('nb_ex_spacers', ls_nb_ex_spacers)
        return ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses

    def get_max_length_loss_sets(self, filter_by=False):
        """
        Returns the maximal length deletions in reconstruction.
        :param filter_by: 'exclude_young_spacers': loss sets are filtered such that only root spacers are contained,
        i.e. only root spacers are considered as deletions for lh ratio computations.
        'contains_old_spacers': Loss sets are filtered such that only sets, that contain at least one old spacer are
        used for lh ratio computations.
        :return:
        """
        if filter_by:
            root_spacers = self.rec_gain_dict[self.root.name]
            if filter_by == 'exclude_young_spacers':
                filter_spacers = set(self.top_order) - set(root_spacers)
            elif filter_by == 'contains_old_spacers':
                filter_spacers = set(root_spacers)
        dict_max_length_losses = {}
        ls_all_max_lengths = []
        ls_all_max_lengths_normalized = []
        for node in self.root.find_clades(order='preorder'):
            possible_sets = copy.deepcopy(node.possible_loss_sets)

            # Could do this after joint sets consolidation (i.e. keeping only joint sets of maximal length)
            if filter_by:
                if filter_by == 'exclude_young_spacers':
                    possible_sets = [p - filter_spacers for p in possible_sets]
                elif filter_by == 'contains_old_spacers':
                    for i, p in enumerate(possible_sets):
                        if not any(f in p for f in filter_spacers):
                            possible_sets.pop(i)
                # failsafe if possible sets is empty to keep standard procedure i.e. no losses: possible_sets = [set()]
                if possible_sets:
                    possible_sets.append(set())

            ls_max_length_losses = []
            ls_max_lengths = []
            while possible_sets:
                lengths_possible_sets = [len(p) for p in possible_sets]
                max_idx = np.argmax(lengths_possible_sets)
                max_length = lengths_possible_sets.pop(max_idx)
                max_set = possible_sets.pop(max_idx)
                if max_length > 0:
                    ls_max_length_losses.append(max_set)
                    ls_max_lengths.append(max_length)
                for i in range(len(possible_sets)):
                    p_set = possible_sets[i]
                    intersection = max_set.intersection(p_set)
                    if intersection:
                        possible_sets[i] = p_set - intersection
            node.max_length_losses = [ls_max_length_losses, ls_max_lengths]
            dict_max_length_losses[node.name] = [ls_max_length_losses, ls_max_lengths]
            if not filter_by:
                self.rec_joint_losses_dict[node.name] = ls_max_length_losses
            ls_all_max_lengths += ls_max_lengths
            if node.up is not None:
                ls_all_max_lengths_normalized += [ml / sum(node.up.rec_spacers) for ml in ls_max_lengths]

        return dict_max_length_losses, ls_all_max_lengths, ls_all_max_lengths_normalized

    def get_rec_duplication_type_counts(self):
        """
        Returns the count of all duplicate candidate event types.
        :return:
        """
        reverse_duplicated_spacers = {}
        for key, value in self.duplicated_spacers.items():
            for v in value:
                reverse_duplicated_spacers[v] = key

        count_duplications = 0
        count_rearrangements = 0
        count_reacquisitions = 0
        count_ind_dups_default = 0
        count_other_dup_events = 0
        count_indep_acquisitions_not_dup_candidates = 0
        set_encountered_gained_spacers = set()

        all_duplicated_spacers = set()
        for value in self.duplicated_spacers.values():
            all_duplicated_spacers.union(value)

        for node in self.root.find_clades():
            # dup_to_spacer_name = [d for d in self.rec_duplications_dict.get(node.name, [])]
            dup_to_unique_name = [reverse_duplicated_spacers[d] for d in self.rec_duplications_dict.get(node.name,
                                                                                                        [])]
            count_duplications += len(dup_to_unique_name)

            # dup_to_spacer_name = [self.spacer_numbers_to_names[d] for d in self.rec_rearrangements_dict.get(node.name,
            #                                                                                                 [])]
            dup_to_unique_name = [reverse_duplicated_spacers[d] for d in self.rec_rearrangements_dict.get(node.name,
                                                                                                          [])]
            count_rearrangements += len(dup_to_unique_name)

            # dup_to_spacer_name = [self.spacer_numbers_to_names[d] for d in ]
            dup_to_unique_name = [reverse_duplicated_spacers[d] for d in self.rec_reacquisition_dict.get(node.name,
                                                                                                         [])]
            count_reacquisitions += len(dup_to_unique_name)

            # dup_to_spacer_name = [self.spacer_numbers_to_names[d]
            #                       for d in ]
            dup_to_unique_name = [reverse_duplicated_spacers[d]
                                  for d in self.rec_indep_gain_dict.get(node.name, [])]
            count_ind_dups_default += len(dup_to_unique_name)

            dup_to_unique_name = [reverse_duplicated_spacers[d] for d in self.rec_other_dup_events_dict.get(node.name,
                                                                                                            [])]
            count_other_dup_events += len(dup_to_unique_name)

            for s in self.rec_gain_dict[node.name]:
                if s in all_duplicated_spacers:
                    continue
                elif s in set_encountered_gained_spacers:
                    count_indep_acquisitions_not_dup_candidates += 1
                else:
                    set_encountered_gained_spacers.add(s)
        return count_duplications, count_rearrangements, count_reacquisitions, count_ind_dups_default, \
            count_other_dup_events, count_indep_acquisitions_not_dup_candidates

    def get_rec_contradictions(self, find_contradictions=True):
        """
        Returns found contradictions. Might be deprecated.
        :param find_contradictions: If True searches for contradictions before returning the dictionary.
        :return:
        """
        if find_contradictions:
            rec_contra_dict = model_tools.find_wrong_order_gains_dict(self.rec_gain_dict, self.root,
                                                                      save_to_tree=True,
                                                                      given_upgraph_order=self.spacer_up_orders)
            self.rec_contra_dict = rec_contra_dict
        return self.rec_contra_dict

    def get_relative_loss_positions(self):
        """
        Returns the relative loss positions, lost spacers, nb of still existing spacers, and dictionaries giving the
        presence of spacers around the first, last (parent) spacer and first equal and last equal spacers.
        :return:
        """
        dict_relative_pos = {}
        leaf_spacers = [leaf.rec_spacers for leaf in self.root.get_terminals()]
        reversed_leaf_spacers = [leaf.rec_spacers[::-1] for leaf in self.root.get_terminals()]

        global_fes = None
        global_les = None
        for i in range(len(leaf_spacers[0])):
            if all(l_s[i] == 1 for l_s in leaf_spacers):
                global_fes = i
                break
        for i in range(len(leaf_spacers[0])):
            if all(l_s[i] == 1 for l_s in reversed_leaf_spacers):
                global_les = i
                break
        global_last_s = len(leaf_spacers[0]) - 1
        for c in self.root.find_clades():
            if c.up is None:
                dict_relative_pos[c.name] = []
                continue
            parent_profile = c.up.rec_spacers
            child_profile = c.rec_spacers
            nb_ex_spacers = sum(parent_profile)
            losses = self.rec_loss_dict[c.name]
            _, loss_idx_pos = differences_child_parent_pos(child_profile, parent_profile)
            loss_pos = [sum(parent_profile[:pos]) for pos in loss_idx_pos]

            global_fes_presence = {}
            m = 0
            if global_fes is not None:
                for i, (p_p, c_p) in enumerate(zip(parent_profile, child_profile)):
                    if i >= global_fes:
                        if p_p == 1:
                            global_fes_presence[m] = True if c_p == 1 else False
                            m += 1

            global_les_presence = {}
            m = 0
            if global_les is not None:
                for i, (p_p, c_p) in enumerate(zip(reversed(parent_profile), reversed(child_profile))):
                    if i >= global_les:
                        if p_p == 1:
                            global_les_presence[m] = True if c_p == 1 else False
                            m += 1

            global_l_s_presence = {}
            m = 0
            for (p_p, c_p) in enumerate(zip(reversed(parent_profile), reversed(child_profile))):
                if p_p == 1:
                    global_l_s_presence[m] = True if c_p == 1 else False
                m += 1

            f_s_presence = {}
            m = 0
            for i, (p_p, c_p) in enumerate(zip(parent_profile, child_profile)):
                if p_p == 1:
                    f_s_presence[m] = True if c_p == 1 else False
                    m += 1

            l_s_presence = {}
            m = 0
            for i, (p_p, c_p) in enumerate(zip(reversed(parent_profile), reversed(child_profile))):
                if p_p == 1:
                    l_s_presence[m] = True if c_p == 1 else False
                    m += 1

            fes_m_presence = {}
            m = 0
            fes_found = False
            for i, (p_p, c_p) in enumerate(zip(parent_profile, child_profile)):
                if fes_found:
                    if p_p == 1:
                        fes_m_presence[m] = True if c_p == 1 else False
                        m += 1
                else:
                    if p_p == c_p:
                        fes_m_presence[m] = True
                        m += 1
                        fes_found = True

            les_m_presence = {}
            m = 0
            les_found = False
            for i, (p_p, c_p) in enumerate(zip(reversed(parent_profile), reversed(child_profile))):
                if les_found:
                    if p_p == 1:
                        les_m_presence[m] = True if c_p == 1 else False
                        m += 1
                else:
                    if p_p == c_p:
                        les_m_presence[m] = True
                        m += 1
                        les_found = True
            dict_relative_pos[c.name] = (
                [1.0 if nb_ex_spacers == 1 else val / (nb_ex_spacers - 1) for val in loss_pos], losses[::-1],
                nb_ex_spacers, (fes_m_presence, les_m_presence, l_s_presence, f_s_presence),
                (global_fes_presence, global_les_presence, global_l_s_presence))
        return dict_relative_pos

    def get_relative_loss_stats(self):
        """
        Based on the relative loss positions from get_relative_loss_positions.
        Returns the relative loss positions, losses, nb of existing spacers, avg. nb of existing spacers, and
        more compact dictionaries about the presence of spacers around the first, last (parent) spacer and
        first/last equal spacers.
        :return:
        """
        dict_relative_pos = self.get_relative_loss_positions()
        ls_rel_loss_pos = []
        ls_losses = []
        ls_nb_ex_spacers = []
        tree_fes_m_presence = {}
        tree_les_m_presence = {}
        tree_l_s_m_presence = {}
        tree_f_s_m_presence = {}

        tree_global_fes_m_presence = {}
        tree_global_les_m_presence = {}
        tree_global_l_s_m_presence = {}
        for c in self.root.find_clades():
            if c.up is None:
                continue
            rel_loss_pos, losses, nb_ex_spacers, fes_les_l_s_m_presence, \
                global_fes_les_l_s_m_presence = dict_relative_pos[c.name]

            for key, val in fes_les_l_s_m_presence[0].items():
                if key in tree_fes_m_presence:
                    tree_fes_m_presence[key][0] += 1 if val else 0
                    tree_fes_m_presence[key][1] += 1
                else:
                    tree_fes_m_presence[key] = [1, 1] if val else [0, 1]
            for key, val in fes_les_l_s_m_presence[1].items():
                if key in tree_les_m_presence:
                    tree_les_m_presence[key][0] += 1 if val else 0
                    tree_les_m_presence[key][1] += 1
                else:
                    tree_les_m_presence[key] = [1, 1] if val else [0, 1]

            for key, val in fes_les_l_s_m_presence[2].items():
                if key in tree_l_s_m_presence:
                    tree_l_s_m_presence[key][0] += 1 if val else 0
                    tree_l_s_m_presence[key][1] += 1
                else:
                    tree_l_s_m_presence[key] = [1, 1] if val else [0, 1]

            for key, val in fes_les_l_s_m_presence[3].items():
                if key in tree_f_s_m_presence:
                    tree_f_s_m_presence[key][0] += 1 if val else 0
                    tree_f_s_m_presence[key][1] += 1
                else:
                    tree_f_s_m_presence[key] = [1, 1] if val else [0, 1]

            for key, val in global_fes_les_l_s_m_presence[0].items():
                if key in tree_global_fes_m_presence:
                    tree_global_fes_m_presence[key][0] += 1 if val else 0
                    tree_global_fes_m_presence[key][1] += 1
                else:
                    tree_global_fes_m_presence[key] = [1, 1] if val else [0, 1]
            for key, val in global_fes_les_l_s_m_presence[1].items():
                if key in tree_global_les_m_presence:
                    tree_global_les_m_presence[key][0] += 1 if val else 0
                    tree_global_les_m_presence[key][1] += 1
                else:
                    tree_global_les_m_presence[key] = [1, 1] if val else [0, 1]

            for key, val in global_fes_les_l_s_m_presence[2].items():
                if key in tree_global_l_s_m_presence:
                    tree_global_l_s_m_presence[key][0] += 1 if val else 0
                    tree_global_l_s_m_presence[key][1] += 1
                else:
                    tree_global_l_s_m_presence[key] = [1, 1] if val else [0, 1]
            ls_rel_loss_pos += rel_loss_pos
            ls_losses += losses
            ls_nb_ex_spacers.append(nb_ex_spacers)
        avg_nb_ex_spacers = np.mean(ls_nb_ex_spacers)
        return ls_rel_loss_pos, ls_losses, ls_nb_ex_spacers, avg_nb_ex_spacers, \
            (tree_fes_m_presence, tree_les_m_presence, tree_l_s_m_presence, tree_f_s_m_presence), \
            (tree_global_fes_m_presence, tree_global_les_m_presence, tree_global_l_s_m_presence)

    ####################################################################################################################
    # Import of data and preprocessing.

    def import_simulation_set_as_reconstruction(self, sim_gain_dict, sim_loss_dict):
        """
        Imports insertion and deletion events from a simulated tree.
        :return:
        """
        self.rec_gain_dict = {}
        self.rec_loss_dict = {}
        for clade in self.root.find_clades():
            self.rec_gain_dict[clade.name] = list(
                map(lambda x: self.spacer_names_to_numbers[x], sim_gain_dict[clade.name]))
            self.rec_loss_dict[clade.name] = list(
                map(lambda x: self.spacer_names_to_numbers[x], sim_loss_dict[clade.name]))

            clade.rec_spacers = clade.bin_spacers
        return

    def import_reverse_aligned_spacers(self, df_reverse_aligned, dict_duplicated_spacers, plot_order=False):
        self.duplicated_spacers = dict_duplicated_spacers
        ls_aligned_arrays = [[value for value in ls if value != '-'] for ls in df_reverse_aligned.values]
        spacer_orders, _, self.order_adj_matrix = model_tools.determine_order_dict(ls_aligned_arrays,
                                                                                   save_path=os.path.join(
                                                                                       self.save_path,
                                                                                       self.model_name + '_order'),
                                                                                   plot_order=plot_order,
                                                                                   )
        cols = [df_reverse_aligned[col] for col in df_reverse_aligned.columns]
        uniques = [list(set(col) - {'-'}) for col in cols]
        top_order = [u[0] for u in uniques]
        return df_reverse_aligned, spacer_orders, top_order

    def prepare_place_aligned_data_on_tree(self, df_aligned, raw_top_order, plot_order=True, only_terminals=True,
                                           cols_to_join=None, reverse_walkthrough=False, dict_duplicated_spacers=None):
        """
        Determines order-dict/lists and places spacer arrays at nodes, note that you can place spacers along the
        complete tree (if we want to use simulations as reconstruction). But in that case, the rec_gain/loss dicts
        need to be replaced.
        :param dict_duplicated_spacers:
        :param reverse_walkthrough: Responsible to have the same labeling for duplicates between forward and reverse
        orientation.
        :param cols_to_join: These are the columns that can be joined, even though they were separated by mafft.
        :param only_terminals:
        :param df_aligned:
        :param raw_top_order:
        :param plot_order:
        :return:
        """
        if dict_duplicated_spacers is None:
            df_aligned, (self.spacer_orders, self.spacer_up_orders), \
                self.top_order = self.rename_duplicates_determine_order_based_on_msa(df_aligned,
                                                                                     raw_top_order,
                                                                                     plot_order=plot_order,
                                                                                     cols_to_join=cols_to_join,
                                                                                     reverse_walkthrough=reverse_walkthrough)
        else:
            df_aligned, (self.spacer_orders, self.spacer_up_orders), \
                self.top_order = self.import_reverse_aligned_spacers(df_aligned,
                                                                     dict_duplicated_spacers,
                                                                     plot_order=plot_order)

        considered_clades = self.root.get_terminals() if only_terminals else self.root.find_clades()
        for clade in considered_clades:
            clade.spacers = [val for val in df_aligned.loc[clade.name] if val != '-']

        df = self.rename_spacers_to_numbers(as_int=True, only_terminals=only_terminals,
                                            reverse_walkthrough=reverse_walkthrough)
        if self.sim_as_rec:
            self.import_simulation_set_as_reconstruction(self.sim_gain_loss_dicts[0], self.sim_gain_loss_dicts[1])
        return df_aligned

    def rename_duplicates_determine_order_based_on_msa(self, df_aligned, full_alignment, plot_order=True,
                                                       separator='', cols_to_join=None, reverse_walkthrough=False):
        """
        Finds duplicated spacers and makes them unique by adding ascii-uppercase separated by "separator".
        (Plotting trees with alignment requires capital
        letters (for some reason)). Then determines the order dictionaries/lists and a topological order.
        :param cols_to_join: {dup: ls_columns}, the columns that are joinable due to being artificially
        separated by mafft.
        :param separator: changes the separation between duplicate signifier and name.
        :param df_aligned: dataframe with aligned spacers
        :param full_alignment: topological order
        :param plot_order: if True, the order is plotted. 'graphviz' performs the plotting with graphviz.
        :return:
        """
        # duplication additions determines the added signs that are used to rename duplicates (they are added separated
        # by separator, dup_add_idx is used, if there are more duplicates than added signs (should not occur normally).
        duplication_additions = list(string.ascii_uppercase)
        # change this if alignment via tree works
        if cols_to_join is None:
            mask = pd.DataFrame(full_alignment).duplicated(keep=False)
            dup_add_idx = {}
        else:
            # Since cols_to_join: {'2': [[*], [*,*]]} is composed of connected components, add_idx must be assigned
            # according to connected component, unless only one connected component exists then they are not traversed.
            # (mask = False). Not all duplications are in cols_to_join, for them the standard is done. dup_add_idx[0] is
            # always this standard counter.
            columns = df_aligned.columns
            mask = pd.DataFrame(full_alignment).duplicated(keep=False).astype(int)
            for i, (s, col) in enumerate(zip(full_alignment, columns)):
                if cols_to_join.get(s, False):
                    if len(cols_to_join[s]) < 2:
                        mask[i] = 0
                    else:
                        for j, ls_cols in enumerate(cols_to_join[s]):
                            if col in ls_cols:
                                mask[i] = j + 1

            dup_add_idx = {s: list(reversed(range(len(cols_to_join.get(s, [])))))
            if reverse_walkthrough else list(range(len(cols_to_join.get(s, []))))
                           for s in set(full_alignment)}

        df_top_order = pd.DataFrame(full_alignment)
        self.duplicated_spacers = dict()
        # Find and rename duplicated spacers.
        # iterator = reversed(list(enumerate(mask))) if reverse_walkthrough else enumerate(mask)
        iterator = enumerate(mask)
        for i, m in iterator:
            if m:
                spacer = df_top_order.iloc[i].values[0]
                if cols_to_join is None:
                    add_str = separator + (dup_add_idx.get(spacer, 0) // len(duplication_additions) + 1) * \
                              duplication_additions[dup_add_idx.get(spacer, 0) % len(duplication_additions)]
                else:
                    add_str = separator + (dup_add_idx[spacer][m - 1] // len(duplication_additions) + 1) * \
                              duplication_additions[dup_add_idx[spacer][m - 1] % len(duplication_additions)]
                df_top_order.iloc[i] = spacer + add_str
                # save duplicated spacer to use later.
                self.duplicated_spacers[spacer] = self.duplicated_spacers.get(spacer, set())
                self.duplicated_spacers[spacer].add(df_top_order.iloc[i].values[0])
                df_aligned.iloc[:, i] = [a if a == '-' else df_top_order.iloc[i].values[0] for a in
                                         df_aligned.iloc[:, i]]
                if cols_to_join is None:
                    dup_add_idx[spacer] = dup_add_idx.get(spacer, 0) + 1
        # Move in between spacers between combined columns.
        current_top_order = df_top_order.T.values.tolist()[0] if not df_top_order.empty else []
        current_df_aligned = df_aligned
        duplicates_in_top_order = [item for item, count in collections.Counter(current_top_order).items()
                                   if count > 1]
        for dup in duplicates_in_top_order:
            dup_col_locs = [i for i, v in enumerate(current_top_order) if v == dup]
            ls_dup_row_locs = []
            set_dup_row_locs = set()
            min_col_loc, max_col_loc = min(dup_col_locs), max(dup_col_locs)
            # Find rows with present spacer in columns with duplicated value
            for col in dup_col_locs:
                dup_row_locs = [i for i, v in enumerate(current_df_aligned.iloc[:, col]) if v == dup]
                ls_dup_row_locs.append(dup_row_locs)
                set_dup_row_locs = set_dup_row_locs.union(set(dup_row_locs))
            max_in_between_spacers = 0
            max_in_between_spacers_col = min_col_loc
            # find spacers in between min and max duplicate columns, find out which has the most in between spacers.
            ls_ibs_cols_left, ls_ibs_cols_right = [], []
            for i, r_locs in enumerate(ls_dup_row_locs):
                left, right = [], []
                set_found_spacers = set()
                for r in r_locs:
                    for j, v in enumerate(current_df_aligned.iloc[r, :]):
                        if v not in {dup, '-'} and v not in set_found_spacers:
                            if min_col_loc <= j <= max_col_loc:
                                if j < dup_col_locs[i]:
                                    left.append(current_df_aligned.iloc[:, j].name)
                                else:
                                    right.append(current_df_aligned.iloc[:, j].name)
                                set_found_spacers.add(v)

                ibs_total_length = len(set_found_spacers)
                if dup_col_locs[i] in {min_col_loc, max_col_loc} and ibs_total_length > max_in_between_spacers:
                    max_in_between_spacers = ibs_total_length
                    max_in_between_spacers_col = dup_col_locs[i]
                ls_ibs_cols_left.append(left)
                ls_ibs_cols_right.append(right)
            # Move all duplicates to max_in_between_spacers_column
            for r in set_dup_row_locs:
                current_df_aligned.iloc[r, max_in_between_spacers_col] = dup
            # Move in between columns that need to be moved
            cols = list(current_df_aligned.columns.values)
            new_cols = []
            new_current_top_order = []

            if max_in_between_spacers_col == min_col_loc:
                # Move all lefts to left of min_col_loc
                set_combined_lefts = set().union(*ls_ibs_cols_left)
                for i, c in enumerate(cols):
                    if c in set_combined_lefts:
                        new_cols.insert(max_in_between_spacers_col, c)
                        new_current_top_order.insert(max_in_between_spacers_col, current_top_order[i])
                    else:
                        new_cols.append(c)
                        new_current_top_order.append(current_top_order[i])

                current_df_aligned = current_df_aligned[new_cols]
                first_found = False
                idx_to_drop = []
                col_to_drop = []
                for i, v in enumerate(new_current_top_order):
                    if v == dup:
                        if not first_found:
                            first_found = True
                        else:
                            idx_to_drop.append(i)
                            col_to_drop.append(current_df_aligned.iloc[:, i].name)
                current_top_order = [v for i, v in enumerate(new_current_top_order) if i not in idx_to_drop]
                current_df_aligned = current_df_aligned.drop(columns=col_to_drop)

            elif max_in_between_spacers_col == max_col_loc:
                set_combined_rights = set().union(*ls_ibs_cols_right)
                current_top_order = current_top_order[::-1]
                max_in_between_spacers_col = len(cols) - max_in_between_spacers_col - 1
                for i, c in enumerate(reversed(cols)):
                    if c in set_combined_rights:
                        new_cols.insert(max_in_between_spacers_col, c)
                        new_current_top_order.insert(max_in_between_spacers_col, current_top_order[i])
                    else:
                        new_cols.append(c)
                        new_current_top_order.append(current_top_order[i])
                current_df_aligned = current_df_aligned[new_cols]
                first_found = False
                idx_to_drop = []
                col_to_drop = []
                for i, v in enumerate(new_current_top_order):
                    if v == dup:
                        if not first_found:
                            first_found = True
                        else:
                            idx_to_drop.append(i)
                            col_to_drop.append(current_df_aligned.iloc[:, i].name)

                current_top_order = [v for i, v in enumerate(new_current_top_order) if i not in idx_to_drop]
                current_df_aligned = current_df_aligned.drop(columns=col_to_drop)
                current_top_order = current_top_order[::-1]
                current_df_aligned = current_df_aligned.iloc[:, ::-1]
            else:
                raise ValueError('Something went wrong...')

        ls_aligned_arrays = [[value for value in ls if value != '-'] for ls in current_df_aligned.values]
        spacer_orders, _, self.order_adj_matrix = model_tools.determine_order_dict(ls_aligned_arrays,
                                                                                   save_path=os.path.join(
                                                                                       self.save_path,
                                                                                       self.model_name + '_order'),
                                                                                   plot_order=plot_order,
                                                                                   )
        return current_df_aligned, spacer_orders, current_top_order

    def rename_spacers_to_numbers(self, order='level', only_terminals=False, as_int=True,
                                  reverse_walkthrough=False):
        """
        Renames the spacers to numbers for internal use. The orders and duplicate candidates are also relabeled
        accordingly.
        """
        # generate enough numbers
        all_spacers_num = list(range(1, len(self.top_order) + 1)) if reverse_walkthrough \
            else list(reversed(range(1, len(self.top_order) + 1)))

        # create relabeling dictionaries
        self.spacer_names_to_numbers = {spacer_name: spacer_nb for spacer_name, spacer_nb in zip(self.top_order,
                                                                                                 all_spacers_num)}
        self.spacer_numbers_to_names = {spacer_nb: spacer_name for spacer_name, spacer_nb in
                                        self.spacer_names_to_numbers.items()}

        # relabel orders
        new_spacer_orders = {}
        new_spacer_up_orders = {}
        for key, order_set in self.spacer_orders.items():
            new_spacer_orders[self.spacer_names_to_numbers[key]] = set(map(lambda x: self.spacer_names_to_numbers[x],
                                                                           order_set))
        for key, order_set in self.spacer_up_orders.items():
            new_spacer_up_orders[self.spacer_names_to_numbers[key]] = set(
                map(lambda x: self.spacer_names_to_numbers[x],
                    order_set))
        self.spacer_orders = new_spacer_orders
        self.spacer_up_orders = new_spacer_up_orders

        # relabel duplicates to numbers.
        new_duplicated_spacers = {}
        for key, duplicates in self.duplicated_spacers.items():
            new_duplicated_spacers[key] = set(map(lambda x: self.spacer_names_to_numbers[x], duplicates))
        self.duplicated_spacers = new_duplicated_spacers

        aligned_spacers = pd.DataFrame(columns=all_spacers_num)

        considered_clades = self.get_terminals(order=order) if only_terminals else self.find_clades(order=order)

        # creates spacer profiles along the tree.
        for c in considered_clades:
            num_spacers = [self.spacer_names_to_numbers[s] for s in c.spacers]

            c.bin_spacers = np.isin(all_spacers_num, test_elements=num_spacers).astype(int) if as_int else np.isin(
                all_spacers_num, test_elements=num_spacers)
            aligned_spacers.loc[c.name] = c.bin_spacers

        aligned_spacers = aligned_spacers.astype(int) if as_int else aligned_spacers
        return aligned_spacers

    ####################################################################################################################
    # Treatment of duplicate candidates.
    def duplication_rearrangement_placement(self, guide_rec_gain_dict):
        """
        Employed after guide reconstruction. Finds the events that are placed closest to the leafs among the candidates
        of each found duplicate. These events are removed from the spacer orders, so they do not produce contradictions
        and are not pulled up by the refinement.
        :return:
        """

        if self.duplicated_spacers is None:
            print('Run msa to detect duplicated spacers first.')
        for spacer_name, set_individual_spacers in self.duplicated_spacers.items():
            encountered_candidates = 0
            for c in self.root.find_clades(order='postorder'):
                if encountered_candidates >= len(set_individual_spacers) - 1:
                    break
                guide_gain_inter = set_individual_spacers.intersection(set(guide_rec_gain_dict.get(c.name, [])))
                while encountered_candidates < len(set_individual_spacers) - 1 and len(guide_gain_inter) > 0:
                    encountered_candidates += 1

                    # dropping
                    spacer_to_remove = guide_gain_inter.pop()
                    self.spacer_orders[spacer_to_remove] = set()
                    self.spacer_up_orders[spacer_to_remove] = set()
                    for spacer in self.top_order:
                        self.spacer_orders[self.spacer_names_to_numbers[spacer]] = self.spacer_orders[
                                                                                       self.spacer_names_to_numbers[
                                                                                           spacer]] - {
                                                                                       spacer_to_remove}
                        self.spacer_up_orders[self.spacer_names_to_numbers[spacer]] = self.spacer_up_orders[
                                                                                          self.spacer_names_to_numbers[
                                                                                              spacer]] - {
                                                                                          spacer_to_remove}

                    self.rec_dup_rearr_candidates[c.name] = self.rec_dup_rearr_candidates.get(c.name, []) \
                                                            + [spacer_to_remove]
        return self.rec_dup_rearr_candidates

    def distinction_dup_rearrangement(self):
        """
        Employed after refinement of reconstruction, distinguishes between rearrangements, duplications, reacquisitions,
        independent acquisitions and others, based on heuristic described in the paper.
        :return:
        """
        for original_name, set_individual_spacers in self.duplicated_spacers.items():
            for node in self.root.find_clades(order='postorder'):
                if node.name in self.rec_dup_rearr_candidates:
                    set_candidates = set(self.rec_dup_rearr_candidates[node.name])
                    candidates = set_candidates.intersection(set_individual_spacers)
                    other_individuals = set_individual_spacers - candidates

                    bool_duplication = False
                    if node.up is None:
                        parent_rec_spacers = np.zeros(len(node.rec_spacers))
                    else:
                        parent_rec_spacers = node.up.rec_spacers
                    parent_existent_spacers, _ = self.differences_child_parent(parent_rec_spacers,
                                                                               np.zeros(len(node.rec_spacers)))
                    existent_spacers, _ = self.differences_child_parent(node.rec_spacers,
                                                                        np.zeros(len(node.rec_spacers)))
                    ls_parent_min_existent_spacer_pos = [i for i, v in enumerate(parent_rec_spacers) if v == 1]
                    parent_min_existent_spacer_pos = min(ls_parent_min_existent_spacer_pos) \
                        if ls_parent_min_existent_spacer_pos else len(parent_rec_spacers)

                    if len(other_individuals.intersection(set(existent_spacers))) > 0:
                        bool_duplication = True

                    path_to_candidate = list(self.root.get_path(node))[:-1]
                    bool_double_gain = False
                    for path_node in path_to_candidate:
                        if len(other_individuals.intersection(set(self.rec_loss_dict[path_node.name]))) > 0:
                            bool_double_gain = True
                    bool_is_acquisition = {}
                    for s in candidates:
                        candidate_pos = len(node.rec_spacers)
                        for idx in range(len(node.rec_spacers)):
                            if len(node.rec_spacers) - idx == s:
                                candidate_pos = idx
                        if candidate_pos < parent_min_existent_spacer_pos:
                            bool_is_acquisition[s] = True
                        else:
                            bool_is_acquisition[s] = False

                    # rearrangement (spacer moved on branch)
                    if len(other_individuals.intersection(set(self.rec_loss_dict[node.name]))) > 0:
                        self.rec_rearrangements_dict[node.name] = self.rec_rearrangements_dict.get(node.name, []) + \
                                                                  list(candidates)
                    # spacer is duplicated, i.e. spacer gain while other exists (on the same branch)
                    elif bool_duplication:
                        for c in candidates:
                            self.rec_duplications_dict[node.name] = self.rec_duplications_dict.get(node.name,
                                                                                                   []) + [c]
                    # spacer gained -> lost -> gained again
                    elif bool_double_gain:
                        for c in candidates:
                            if bool_is_acquisition.get(c, False):
                                self.rec_reacquisition_dict[node.name] = self.rec_reacquisition_dict.get(node.name, []) \
                                                                         + [c]
                            else:
                                # Default case
                                self.rec_other_dup_events_dict[
                                    node.name] = self.rec_other_dup_events_dict.get(node.name, []) + [c]
                    # this is the default case (most likely an independent acquisition)
                    else:
                        for c in candidates:
                            if bool_is_acquisition.get(c, False):
                                self.rec_indep_gain_dict[node.name] = self.rec_indep_gain_dict.get(node.name, []) + [c]
                            else:
                                self.rec_other_dup_events_dict[node.name] = self.rec_other_dup_events_dict.get(
                                    node.name,
                                    []) + [c]
        return self.rec_duplications_dict, self.rec_rearrangements_dict, self.rec_reacquisition_dict, \
            self.rec_indep_gain_dict, self.rec_other_dup_events_dict

    ####################################################################################################################
    # generate data for lh ratio test, and compute lh ratio.
    def generate_possible_losses_sets(self):
        dict_matrices = {}
        for node in self.root.find_clades(order='preorder'):
            if node.up is None:
                node.possible_loss_sets = [set()]
                dict_matrices[node.name] = node.possible_loss_sets
                continue
            node.possible_loss_sets = self.generate_possible_loss_set(node)
            dict_matrices[node.name] = node.possible_loss_sets
        return dict_matrices

    def generate_possible_loss_set(self, node):
        if node.up is None:
            raise ValueError('Given node is root!')

        parent_profile = node.up.rec_spacers
        child_profile = node.rec_spacers
        _, losses = self.differences_child_parent(child_profile, parent_profile)
        a = [set(losses)]  # or all orders?
        l, r = self.get_l_r_sets_order_dict(node)
        # print(node.name, 'l', l, 'r', r)
        for l_k, r_k in zip(l, r):
            new_a = []
            for a_i in a:
                a_i_l = a_i.intersection(l_k)
                a_i_r = a_i.intersection(r_k)
                a_i_v = a_i - a_i_l.union(a_i_r)

                if len(a_i_l) > 0 and len(a_i_r) > 0:
                    new_a.append(a_i_l.union(a_i_v))
                    new_a.append(a_i_r.union(a_i_v))
                else:
                    new_a.append(a_i)
            a = new_a
        return a

    def get_l_r_sets_order_dict(self, node):
        if node.up is None:
            raise ValueError('Given node is root!')
        parent_profile = node.up.rec_spacers
        child_profile = node.rec_spacers
        _, losses = self.differences_child_parent(child_profile, parent_profile)
        losses = set(losses)

        ex_spacers = prof2seq(
            np.array([1 if (p == 1 and c == 1) else 0 for c, p in zip(child_profile, parent_profile)]))
        l, r = [], []
        for k in ex_spacers:
            l.append(self.spacer_up_orders[k].intersection(losses))
            r.append(self.spacer_orders[k].intersection(losses))
        return l, r

    def compute_lh_ratio_for_given_lh_fct(self, method='Nelder-Mead', filter_by=False):
        """
        Function to optimize lh ratio of ilm vs blm, note that ilm likelihood is computed with old lh fct (should make
        no difference, but new lh function can produce artifacts).
        :param method:
        :param filter_by:
        :return:
        """
        ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses = self.get_data_for_lh_ratio(filter_by=filter_by)
        lh_fct_0 = lambda x: -compute_tree_lh(x, 1.0, ls_bl, ls_nb_ex_spacers,
                                              ls_nb_max_length_losses,
                                              log=True)
        lh_fct_1 = lambda x: -self.compute_tree_lh_for_given_lh_fct(x[0], x[1], ls_bl, ls_nb_ex_spacers,
                                                                    ls_nb_max_length_losses,
                                                                    log=True)
        # start_values_0 = 0.5
        start_values_1 = np.array([0.5, 1.5])
        bnds_0 = (0.0, 100000)
        bnds_1 = ((0.0, None), (1.0, None))

        max_0 = minimize_scalar_fct(lh_fct_0, bounds=bnds_0, method='bounded')
        max_1 = minimize_fct(lh_fct_1, start_values_1, bounds=bnds_1, method=method)

        lh_0 = np.exp(-max_0.fun)
        lh_1 = np.exp(-max_1.fun)

        return lh_0, lh_1, max_0, max_1

    def compute_tree_lh_for_given_lh_fct(self, loss_rate, alpha, ls_bl, nb_survivors, nb_max_length_losses, log=True):
        """
        Computes bdm-lh of a tree using a given likelihood function.
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
            max_mll = max(mll) if mll else 0

            ls_f = np.log([max(f(bl, loss_rate, alpha), self.eps) for f in self.lambdifyed_lh_fct[:max_mll + 1]])
            ls_f[np.isnan(ls_f) | np.isinf(ls_f)] = -self.big_eps
            ln_keep = - loss_rate * bl * n

            # Handle blocks that are longer than the available likelihood functions. Divide them in max length blocks
            # + Rest.
            if max_mll + 1 > len(ls_f):
                ls_ln_mult_loss_lh = []
                for k_i in mll:
                    val = 0
                    val += ls_f[k_i % (len(ls_f) - 1)]
                    incl_number = k_i // (len(ls_f) - 1)
                    for j in range(incl_number):
                        val = incl_number * ls_f[-1]
                    ls_ln_mult_loss_lh.append(val)
            else:
                ls_ln_mult_loss_lh = [ls_f[k_i] for k_i in mll]
            ln_mult_loss_lh = sum(ls_ln_mult_loss_lh)
            tree_lh.append(ln_mult_loss_lh + ln_keep)
        return sum(tree_lh) if log else np.exp(sum(tree_lh))

    def compute_lh_ratio(self, method='Nelder-Mead', filter_by=False, alpha_bias_correction=False,
                         rho_bias_correction=False):
        """
        Computes the "naive" lh ratio of the tree with optional bias corrections.
        :param method:
        :param filter_by:
        :param alpha_bias_correction:
        :param rho_bias_correction:
        :return:
        """
        ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses = self.get_data_for_lh_ratio(filter_by=filter_by)
        lh_fct_0 = lambda x: -compute_tree_lh(x, 1, ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses,
                                              log=True,
                                              rec_tree=None, gains=None)
        lh_fct_1 = lambda x: -compute_tree_lh(x[0], x[1], ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses, log=True,
                                              rec_tree=None, gains=None)
        # start_values_0 = 0.5
        start_values_1 = np.array([0.5, 1.5])
        bnds_0 = (0, 100000)
        bnds_1 = ((0, None), (1, None))

        max_0 = minimize_scalar_fct(lh_fct_0, bounds=bnds_0, method='bounded')
        max_1 = minimize_fct(lh_fct_1, start_values_1, bounds=bnds_1, method=method)

        if rho_bias_correction:
            tree = self
            gains = self.gains_as_ls_for_lh_computation()
            lh_fct_bias_corr = lambda x: -compute_tree_lh(x, 1, ls_bl, ls_nb_ex_spacers, ls_nb_max_length_losses,
                                                          log=True,
                                                          rec_tree=tree, gains=gains)
            max_0_bias_corr = minimize_scalar_fct(lh_fct_bias_corr, bounds=bnds_0, method='bounded')
            rho_correction = max_0_bias_corr.x - max_0.x
            max_1.x[0] = max_1.x[0] + rho_correction

            max_0.x = max_0_bias_corr.x
        if alpha_bias_correction:
            max_1.x[1] = self.alpha_bias_correction(max_1.x[1], ls_nb_max_length_losses)
        lh_0 = np.exp(-max_0.fun)
        lh_1 = np.exp(-max_1.fun)

        return lh_0, lh_1, max_0, max_1

    def compute_unobserved_lh(self, eff_lr):
        """
        Computes the likelihoods (along the tree) of a gain to be unobserved in all leafs according to the IDM.
        :param eff_lr:
        :return:
        """
        ls_unobserved_lh = []
        for c in self.root.find_clades(order='postorder'):
            if c.up is None:
                c.unobserved_lh = np.sum([d.unobserved_lh for d in c.clades])
            elif c.is_terminal():
                c.unobserved_lh = 1 - np.exp(-eff_lr * c.branch_length)
                c.unobserved_lh = np.log(c.unobserved_lh) if c.unobserved_lh < 1 else np.log(1 - 1e-16)
            else:
                c.unobserved_lh = 1 - np.exp(-eff_lr * c.branch_length)
                c.unobserved_lh += np.exp(-eff_lr * c.branch_length + np.sum([d.unobserved_lh for d in c.clades]))
                c.unobserved_lh = np.log(c.unobserved_lh) if c.unobserved_lh < 1 else np.log(1 - 1e-16)
        for c in self.root.find_clades(order='preorder'):
            if c.is_terminal():
                ls_unobserved_lh.append(None)
            else:
                ls_unobserved_lh.append(np.sum([d.unobserved_lh for d in c.clades]))  # c.unobserved_lh
        return ls_unobserved_lh

    def gains_as_ls_for_lh_computation(self):
        """
        Returns the number of gains for each branch in a list.
        :return:
        """
        ls_gains = []
        for c in self.root.find_clades(order='preorder'):
            if c.is_terminal():
                ls_gains.append(0)
            ls_gains.append(len(self.rec_gain_dict[c.name]))
        return ls_gains

    def alpha_bias_correction(self, alpha_est, ls_nb_mll):
        """
        Computes the alpha bias correction for a given list of max. length losses.
        :param ls_nb_mll:
        :return:
        """
        p_est = 1 / alpha_est
        nb_events = sum([len(nb_mll) for nb_mll in ls_nb_mll])
        sqrt = np.sqrt((nb_events + 1) ** 2 - 4 * nb_events * p_est)
        bias = -0.5 * (sqrt - nb_events + 2 * p_est - 1)
        return max(1, 1 / (p_est + bias))

    ####################################################################################################################
    # Reconstruction functions.
    def ml_anc_joint_correct_contra(self, repeated=False, finite_bl_at_root=True):
        """
        Returns logarithm (!) of the root joint max. likelihood and of the max. likelihood of each state.
        :param finite_bl_at_root: If True, runs algorithm where the maximum of clade branch lengths is used for the
        model at root (instead of an infinite branch length (stationary model)).
        :param repeated: :return:
        """
        joint_lh, seq_lh = self.ml_anc_joint(fixed_states=False)
        contra_dict = self.get_rec_contradictions(find_contradictions=True)
        for clade in self.root.find_clades():
            clade.fixed_cx = []
        while contra_dict:
            for clade in self.root.find_clades(order='postorder'):
                if clade.up is None:
                    continue
                for contra in clade.contra_gain:
                    contra_pos = len(clade.rec_spacers) - contra
                    fixed_cx = np.array([1, 1])
                    clade.up.fixed_cx.append([contra_pos, fixed_cx])
            joint_lh, seq_lh = self.ml_anc_joint(fixed_states=True, finite_bl_at_root=finite_bl_at_root)
            contra_dict = self.get_rec_contradictions(find_contradictions=True)
            if not repeated:
                contra_dict = {}
        return joint_lh, seq_lh

    def ml_anc_joint(self, fixed_states=False, finite_bl_at_root=True):
        """
        Adapted from TreeTime https://github.com/neherlab/treetime.
        :param finite_bl_at_root: If False use stationary distribution, otherwise use max of branch lengths of
        subclades and use standard transition probabilities.
        :param fixed_states:
        :return:
        """
        L = len(self.top_order)
        n_states = 2

        for node in self.root.find_clades(order='postorder'):
            if node.up is None:
                node.joint_cx = None
                continue
            branch_length = node.branch_length
            log_transitions = np.log(np.maximum(self.eps, self.rec_model.prob(branch_length))).T
            if node.is_terminal():
                tmp_prof = seq2prof(node.bin_spacers)
                msg_from_children = np.log(np.maximum(tmp_prof, self.eps))
                msg_from_children[np.isnan(msg_from_children) | np.isinf(msg_from_children)] = -self.big_eps
            else:
                msg_from_children = np.sum(np.stack([c.joint_lx for c in node.clades], axis=0), axis=0)

            node.joint_lx = np.zeros((L, n_states))
            node.joint_cx = np.zeros((L, n_states), dtype=int)
            for char_idx, char in enumerate(range(2)):
                # Pij(idx)*L_ch(idx)
                msg_to_parent = (log_transitions[:, char_idx].T + msg_from_children)
                # which char realizes the maximum
                node.joint_cx[:, char_idx] = msg_to_parent.argmax(axis=1)
                # actual maximum likelihood, axis=1: max of each row
                node.joint_lx[:, char_idx] = msg_to_parent.max(axis=1)
                if fixed_states:
                    if node.fixed_cx:
                        for pos, state in node.fixed_cx:
                            node.joint_cx[pos, char_idx] = state[char_idx]
                            node.joint_lx[pos, char_idx] = msg_to_parent[pos, state[char_idx]]

        msg_from_children = np.sum(np.stack([c.joint_lx for c in self.root], axis=0), axis=0)
        # Pi (stationary distribution) * prod_ch l_ch(i)
        # contribution of stationary distribution
        if finite_bl_at_root:
            max_bl = max([clade.branch_length for clade in self.root.clades])
            log_transitions = np.log(np.maximum(self.eps, self.rec_model.prob(max_bl)[0, :]))
            self.root.joint_lx = msg_from_children + log_transitions
        else:
            self.root.joint_lx = msg_from_children + np.log(self.rec_model.stationary_p[0, :])

        normalized_profile = (self.root.joint_lx.T - self.root.joint_lx.max(axis=1)).T

        # choose sequence characters from profile, different from general alphabets
        seq = prof2seq(np.exp(normalized_profile))
        if fixed_states:
            if self.root.fixed_cx:
                for pos, state in self.root.fixed_cx:
                    seq[pos] = state[0]

        self.root.sequence_lh = np.choose(seq, self.root.joint_lx.T)
        self.root.sequence_joint_lh = self.root.sequence_lh.sum()
        self.root.rec_spacers = seq
        len_spacers = len(seq)
        self.rec_gain_dict[self.root.name], self.rec_loss_dict[self.root.name] = self.differences_child_parent(
            self.root.rec_spacers, np.zeros(len_spacers))

        nodes_to_reconstruct = self.root.find_clades(order='preorder')
        for node in nodes_to_reconstruct:
            if node.up is None:
                continue
            node.rec_spacers = np.choose(node.up.rec_spacers, node.joint_cx.T)
            self.rec_gain_dict[node.name], self.rec_loss_dict[node.name] = self.differences_child_parent(
                node.rec_spacers,
                node.up.rec_spacers)
        return self.root.sequence_joint_lh, self.root.sequence_lh

    ####################################################################################################################
    # Visualization functions.
    def visualize_order(self, save_folder, graph_name, save_path_dot, color_dict=None, spacer_or_number='number',
                        provided_numbering=None):
        adj_matrix, label_dict = self.order_adj_matrix
        if spacer_or_number == 'number':
            if provided_numbering is None:
                label_dict = {k: str(self.spacer_names_to_numbers[str(v)]) for k, v in label_dict.items()}
            else:
                label_dict = {k: str(provided_numbering[str(v)]) for k, v in label_dict.items()}

        vis.plot_order_w_graphviz(adj_matrix, label_dict=label_dict, do_show=False,
                                  save_folder=save_folder,
                                  graph_name=graph_name,
                                  file_name=save_path_dot,
                                  color_dict=color_dict)

    def visualize_tree(self, name='default', do_show=True, path='pictures', sort_top_order='nb_leafs',
                       determine_fsizes_by_str_len=True, indicate_joint_del=False, provided_numbering=None,
                       provided_bg_colors=None, spacer_labels_num=True, re_plot_order=True):
        df = pd.DataFrame(columns=['rec_spacers'])
        df_changes = pd.DataFrame(columns=['rec_gains', 'rec_losses'])
        rec_losses_dict = self.rec_joint_losses_dict if indicate_joint_del else self.rec_loss_dict

        top_order_sorted_by_age = []
        ls_nb_leafs = []
        ls_nodes = []
        # get node iterator
        if isinstance(sort_top_order, str):
            if sort_top_order == 'nb_leafs':
                for node in self.root.find_clades():
                    nb_leafs = len(node.get_terminals())
                    ls_nb_leafs.append(nb_leafs)
                    ls_nodes.append(node)
                ls_sorted_nb_leafs, ls_sorted_nodes = (list(t) for t in zip(*sorted(zip(ls_nb_leafs, ls_nodes),
                                                                                    key=lambda x: x[0])))
                ls_sorted_nb_leafs, ls_sorted_nodes = list(reversed(ls_sorted_nb_leafs)), list(
                    reversed(ls_sorted_nodes))
                node_iterator = ls_sorted_nodes  # self.root.find_clades(order=order)
            elif sort_top_order in {'preorder', 'postorder', 'level'}:
                node_iterator = self.root.find_clades(order=sort_top_order)
            else:
                node_iterator = self.root.find_clades()
        else:
            node_iterator = self.root.find_clades()

        # sort by iterator order
        for node in node_iterator:
            branch_gains = [self.spacer_numbers_to_names[s] for s in self.rec_gain_dict[node.name]]
            top_order_sorted_by_age = list(reversed(branch_gains)) + top_order_sorted_by_age

            # handling duplication candidates. Probably better to handle this by using order, but candidates are removed
            # from order... Would be a bit messy. Should at least respect order within leaf arrays now!
            for s in self.rec_dup_rearr_candidates.get(node.name, []):
                s_as_name = self.spacer_numbers_to_names[s]
                idx_in_top_order_sorted = [i for i, v in enumerate(top_order_sorted_by_age)
                                           if v == s_as_name][0]
                spacers_before = set()
                for leaf in self.root.get_terminals():
                    if s_as_name in leaf.spacers:
                        idx = [i for i, v in enumerate(leaf.spacers) if v == s_as_name][-1]
                        spacers_before = spacers_before.union(set(leaf.spacers[:idx]))

                spacers_behind_dup_cand = set(top_order_sorted_by_age[idx_in_top_order_sorted + 1:])
                intersection = spacers_behind_dup_cand.intersection(spacers_before)
                if intersection:
                    insert_loc = max([i for i, v in enumerate(top_order_sorted_by_age) if v in intersection])
                    spacer = top_order_sorted_by_age.pop(idx_in_top_order_sorted)
                    top_order_sorted_by_age.insert(insert_loc, spacer)
            df_changes.loc[node.name] = [self.rec_gain_dict[node.name],
                                         rec_losses_dict[node.name]]
        top_order_sorted_by_age = [s for s in self.top_order if
                                   s not in top_order_sorted_by_age] + top_order_sorted_by_age

        # If multiple gains happen uniqify the topological order.
        if len(top_order_sorted_by_age) > len(self.top_order):
            top_order_sorted_by_age = list(
                reversed(misc.uniqify_ls_order_preserving(top_order_sorted_by_age[::-1])))

        if not isinstance(sort_top_order, str):
            top_order_sorted_by_age = sort_top_order

        if sort_top_order == 'top_order':
            top_order_sorted_by_age = self.top_order
        name_top_order = {i: s for i, s in enumerate(self.top_order)}
        pos_top_order_sorted_by_age = {s: i for i, s in enumerate(top_order_sorted_by_age)}

        # resort reconstructed spacers
        for node in self.root.find_clades(order='level'):
            new_array = list(np.zeros(len(node.rec_spacers)).astype(int))
            for i, s in enumerate(node.rec_spacers):
                if s == 1:
                    spacer_name = name_top_order[i]
                    spacer_pos = pos_top_order_sorted_by_age[spacer_name]
                    new_array[spacer_pos] = 1
            df.loc[node.name] = [new_array]

        # self.write_vis_data_to_json(df, df_changes,
        #                             {'rec_contra_dict': self.rec_contra_dict,
        #                              'rec_duplications_dict': self.rec_duplications_dict,
        #                              'rec_rearrangements_dict': self.rec_rearrangements_dict,
        #                              'rec_reacquisition_dict': self.rec_reacquisition_dict,
        #                              'rec_indep_gain_dict': self.rec_indep_gain_dict,
        #                              'rec_other_dup_events_dict': self.rec_other_dup_events_dict, },
        #                             self.spacer_names_to_numbers,
        #                             top_order_sorted_by_age,
        #                             name, path)
        if re_plot_order:
            self.visualize_order(path, name, save_path_dot=os.path.join(name + '_order.dot'),
                                 color_dict=provided_bg_colors,
                                 spacer_or_number='number' if spacer_labels_num else 'spacer')

        dict_bg_colors = vis.visualize_tree_rec_spacers(str(self.format(fmt='newick')),
                                                        df_rec_spacers=df,
                                                        df_gains_losses=df_changes,
                                                        name=name,
                                                        do_show=do_show,
                                                        rec_contra_dict=self.rec_contra_dict,
                                                        path=path,
                                                        top_order=top_order_sorted_by_age,
                                                        alphabet=self.spacer_names_to_numbers,
                                                        rec_duplications_dict=self.rec_duplications_dict,
                                                        rec_rearrangements_dict=self.rec_rearrangements_dict,
                                                        rec_double_gains_dict=self.rec_reacquisition_dict,
                                                        rec_indep_gains_dict=self.rec_indep_gain_dict,
                                                        rec_other_dup_events_dict=self.rec_other_dup_events_dict,
                                                        determine_fsizes_by_str_len=determine_fsizes_by_str_len,
                                                        indicate_joint_del=indicate_joint_del,
                                                        provided_numbering=provided_numbering,
                                                        provided_bg_colors=provided_bg_colors,
                                                        spacer_labels_num=spacer_labels_num,
                                                        )
        return dict_bg_colors

    def write_vis_data_to_json(self, df, df_changes, dict_other_events, spacer_names_to_numbers,
                               top_order, name, path):
        print('path', path)
        if not os.path.exists(os.path.join(path, name)):
            os.makedirs(os.path.join(path, name))
        with open(os.path.join(path, name, name + '_tree.nwk'), 'w') as f:
            f.write(str(self.format(fmt='newick')))
        df.to_json(os.path.join(path, name, name + '_rec_spacers.json'))
        df_changes.to_json(os.path.join(path, name, name + '_rec_gains_losses.json'))
        with open(os.path.join(path, name, name + '_other_events.json'), 'w') as f:
            json.dump(dict_other_events, f)
        with open(os.path.join(path, name, name + '_spacer_names_to_numbers.json'), 'w') as f:
            json.dump(spacer_names_to_numbers, f)
        with open(os.path.join(path, name, name + '_top_order.json'), 'w') as f:
            json.dump(top_order, f)

    ####################################################################################################################
    # Helper functions.
    def differences_child_parent(self, child_profile, parent_profile):
        diff_p_ch = child_profile - parent_profile
        gains = [self.spacer_names_to_numbers[self.top_order[idx]] for idx, val in enumerate(diff_p_ch) if val == 1]
        losses = [self.spacer_names_to_numbers[self.top_order[idx]] for idx, val in enumerate(diff_p_ch) if val == -1]
        gains = gains[::-1]
        losses = losses[::-1]
        return gains, losses
