import numpy as np
from Bio import Phylo
import copy
import pickle
import itertools
import pandas as pd

from model.data_classes.advanced_tree import AdvancedTree
from model.spacer_visualization import visualization as vis
from model.data_classes.crisprdata import CRISPRGroup, CRISPRArray, Species
from model.helpers import misc, stats
from model import model_tools


class SimulationTree(AdvancedTree):
    """
    Class to simulate a group of CRISPR arrays along a given tree.
    """

    def __init__(self, sim_deletion_model, event_sim, gtr_gain_rate, gtr_loss_rate, deletion_model_parameter,
                 save_path, root, rooted, model_name=None, blm_rate_correction_last_spacer=True,
                 blm_two_directional=False, sim_dup_events=None, *args, **kwargs):
        """

        :param simulation_model_config: available loss models: 'block', 'fragment', 'independent'
        :param save_path:
        :param root:
        :param rooted:
        :param model_name:
        :param args:
        :param kwargs:
        """
        super(SimulationTree, self).__init__(root=root,
                                             rooted=rooted,
                                             model_name=model_name,
                                             *args,
                                             **kwargs)

        self.gain_dict = {}
        self.loss_dict = {}
        self.block_loss_dict = {}
        self.gain_loss_events_dict = {}

        self.eps = 1e-12
        self.big_eps = 1e10

        self.next_spacer = 1

        self.sim_deletion_model = sim_deletion_model
        self.event_sim = event_sim
        self.gain_rate = gtr_gain_rate
        self.loss_rate = gtr_loss_rate

        self.alpha = deletion_model_parameter

        self.save_path = save_path

        self.blm_rate_correction_last_spacer = blm_rate_correction_last_spacer
        self.blm_two_directional = blm_two_directional

        self.sim_dup_events = sim_dup_events

    def get_array_lengths(self):
        dict_lengths = {}
        for leaf in self.root.get_terminals():
            dict_lengths[leaf.name] = len([s for s in leaf.spacers if s != 0])
        return dict_lengths

    def avg_length_leaf_spacers(self):
        ls_spacers_len = []
        for c in self.root.get_terminals():
            ls_spacers_len.append(len([s for s in c.spacers if s != 0]))
        return np.mean(ls_spacers_len)

    def get_never_seen_spacers(self):
        top_order = set(self.get_topological_order())
        set_seen_spacers = set()
        for n in self.root.get_terminals():
            set_seen_spacers = set_seen_spacers.union(set(n.spacers))
        set_never_seen_spacers = top_order - set_seen_spacers
        return set_never_seen_spacers, len(set_never_seen_spacers)

    def get_gain_loss_dicts(self):
        for key, gains in self.gain_dict.items():
            self.gain_dict[key] = [str(s) for s in gains]
        for key, losses in self.loss_dict.items():
            self.loss_dict[key] = [str(s) for s in losses]
        return self.gain_dict, self.loss_dict

    def get_gain_loss_block_event_dicts(self):
        for key, ls_events in self.gain_loss_events_dict.items():
            new_ls_events = []
            for event in ls_events:
                new_ls_events.append(
                    (event[0], event[1], [str(s) for s in event[2]], event[3], event[4], [str(s) for s in event[5]]))
            self.gain_loss_events_dict[key] = new_ls_events
        return self.gain_dict, self.loss_dict, self.block_loss_dict, self.gain_loss_events_dict

    def get_spacer_arrays(self, only_terminals=False):
        considered_clades = self.root.get_terminals() if only_terminals else self.root.find_clades()
        dict_spacer_arrays = {}
        for c in considered_clades:
            dict_spacer_arrays[c.name] = [str(s) for s in c.spacers if s != 0]
        return dict_spacer_arrays

    def get_topological_order(self):
        return list(reversed(range(1, self.next_spacer)))

    def return_simulation_results(self, only_terminals=True):
        dict_spacer_arrays = self.get_spacer_arrays(only_terminals=only_terminals)
        gain_dict, loss_dict, block_loss_dict, gain_loss_event_dict = self.get_gain_loss_block_event_dicts()
        top_order = self.get_topological_order()
        return dict_spacer_arrays, gain_dict, loss_dict, block_loss_dict, gain_loss_event_dict, top_order

    def get_sim_rel_loss_pos(self):
        ls_nb_ex_spacers = []
        ls_times = []
        ls_losses = []
        ls_loss_pos = []
        ls_rel_loss_pos = []
        for c in self.root.find_clades():
            if c.up is None:
                continue
            for event_tuple in c.gain_loss_events:
                if event_tuple[0] == '-':
                    t = event_tuple[1]
                    losses = event_tuple[2]
                    nb_ex_spacers = event_tuple[3]
                    loss_pos = event_tuple[4]
                    rel_loss_pos = [1.0 if nb_ex_spacers == 1 else lp / (nb_ex_spacers - 1) for lp in loss_pos]

                    ls_times.append(t)
                    ls_losses.append(losses)
                    ls_nb_ex_spacers.append(nb_ex_spacers)
                    ls_loss_pos += loss_pos
                    ls_rel_loss_pos += rel_loss_pos
        return ls_times, ls_losses, ls_nb_ex_spacers, ls_loss_pos, ls_rel_loss_pos

    def get_relative_loss_positions_branch_based(self):
        ls_relative_loss_pos = []
        ls_losses = []
        ls_nb_ex_spacers = []

        for c in self.root.find_clades():
            if c.up is None:
                continue
            nb_ex_spacers = len([c for c in c.up.spacers if c != 0])
            losses = self.loss_dict[c.name]
            loss_pos = [i for i, s in enumerate(c.up.spacers) if s in losses]
            # loss_pos = [sum(parent_profile[:pos]) for pos in loss_idx_pos]

            ls_relative_loss_pos.append([1.0 if nb_ex_spacers == 1 else val / (nb_ex_spacers - 1) for val in loss_pos])
            ls_losses.append(losses)
            ls_nb_ex_spacers.append(nb_ex_spacers)
        return ls_relative_loss_pos, ls_losses, ls_nb_ex_spacers

    def ls_spacer_arrays_to_profile(self, ls_arrays, top_order=None):
        top_order = self.get_topological_order() if top_order is None else top_order
        ls_profiles = []
        for array in ls_arrays:
            new_array = [1 if s in array else 0 for s in top_order]
            ls_profiles.append(new_array)
        return ls_profiles

    def get_event_based_fes_les_pos(self):
        return

    def get_branch_based_fes_les_pos(self):
        dict_fes_les = {}
        top_order = self.get_topological_order()
        leaf_spacers = [leaf.spacers for leaf in self.root.get_terminals()]
        reversed_leaf_spacers = [leaf.spacers[::-1] for leaf in self.root.get_terminals()]

        leaf_spacers = self.ls_spacer_arrays_to_profile(leaf_spacers, top_order=top_order)
        reversed_leaf_spacers = self.ls_spacer_arrays_to_profile(reversed_leaf_spacers, top_order=top_order)

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
                dict_fes_les[c.name] = []
                continue
            parent_profile, child_profile = self.ls_spacer_arrays_to_profile([c.up.spacers, c.spacers],
                                                                             top_order=top_order)

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

            l_s_presence = {}
            m = 0
            for i, (p_p, c_p) in enumerate(zip(reversed(parent_profile), reversed(child_profile))):
                if p_p == 1:
                    l_s_presence[m] = True if c_p == 1 else False
                    m += 1

            f_s_presence = {}
            m = 0
            for i, (p_p, c_p) in enumerate(zip(parent_profile, child_profile)):
                if p_p == 1:
                    f_s_presence[m] = True if c_p == 1 else False
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
            dict_fes_les[c.name] = ((fes_m_presence, les_m_presence, l_s_presence, f_s_presence),
                                    (global_fes_presence, global_les_presence, global_l_s_presence))
        return dict_fes_les

    def get_relative_loss_stats(self):
        dict_relative_pos = self.get_branch_based_fes_les_pos()
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
            fes_les_l_s_m_presence, global_fes_les_l_s_m_presence = dict_relative_pos[c.name]

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
        return (tree_fes_m_presence, tree_les_m_presence, tree_l_s_m_presence, tree_f_s_m_presence), \
            (tree_global_fes_m_presence, tree_global_les_m_presence, tree_global_l_s_m_presence)

    def get_data_as_crispr_group(self, name=None, only_terminals=True):
        dict_spacer_arrays = self.get_spacer_arrays(only_terminals=only_terminals)
        ls_crispr_arrays = [CRISPRArray(name, 'sim', 1, 4, 'sim', 'sim',
                                        'sim', 'sim',
                                        array, 'sim', cas_genes=[], all_cas_types=['sim'],
                                        all_cas_genes=None,
                                        species=Species('sim', 'sim', 'sim', 'sim',
                                                        'sim', 'sim', 'sim', 'sim'),
                                        array_start=None,
                                        array_end=None) for name, array in dict_spacer_arrays.items()]
        crispr_group = CRISPRGroup('sim', ls_crispr_arrays, name=name)
        gain_dict, loss_dict = self.get_gain_loss_dicts()
        crispr_group.gain_loss_dicts = (gain_dict, loss_dict)
        crispr_group.top_order = self.get_topological_order()
        return crispr_group

    def save_to_pkl(self):
        with open(self.save_path + '.pkl', 'wb') as f:
            pickle.dump(self, f)
        return

    def sample_stationary_length_model(self):
        if self.sim_deletion_model == 'block':
            # print(self.gain_rate, self.loss_rate, self.alpha)
            return np.random.poisson(self.gain_rate / (self.loss_rate * self.alpha))
        elif self.sim_deletion_model == 'fragment':
            return stats.sample_stationary_frag_length(1, self.gain_rate / self.loss_rate, dist_limit=1000)[0]
        return np.random.poisson(self.gain_rate / self.loss_rate)

    def generate_new_spacers(self, clade_name, nb_gains):
        spacers = list(reversed(range(self.next_spacer, self.next_spacer + nb_gains)))
        self.gain_dict[clade_name] = self.gain_dict.get(clade_name, []) + list(reversed(spacers))
        self.next_spacer += nb_gains
        return spacers

    def sample_gain_events(self, branch_length):
        lam = self.gain_rate / self.loss_rate * (1 - np.exp(-self.loss_rate * branch_length))
        return np.random.poisson(lam=lam)

    def sample_loss_events_ilm(self, c):
        branch_length = c.branch_length
        parent_spacers = c.up.spacers
        n = len(parent_spacers)
        prob = np.exp(-self.loss_rate * branch_length)
        # prob = np.exp(-np.log(self.loss_rate + 1 / self.loss_rate) * branch_length)

        loss_mask = [1 if r < prob else 0 for r in np.random.rand(n)]
        return loss_mask

    def sample_loss_events_blm(self, c):
        branch_length = c.branch_length
        parent_spacers = c.up.spacers
        n = len(parent_spacers)
        prob = np.exp(-self.loss_rate * branch_length)
        # prob = np.exp(-np.log(self.loss_rate + 1 / self.loss_rate) * branch_length)
        p = 1 / self.alpha
        loss_start_prob = 1 - prob
        loss_mask = list(np.ones(n))
        i = 0
        while i < n:
            if parent_spacers[i] == 0:
                i += 1
                continue
            non_zero_parent_spacers = [True if j != 0 else False for j in parent_spacers[i:]]
            nb_non_zero_parent_spacers = sum(non_zero_parent_spacers)
            uni = np.random.uniform(low=0.0, high=1.0)
            if uni < loss_start_prob:
                fragment_length = np.random.geometric(p)
                fragment_length = min(fragment_length, nb_non_zero_parent_spacers)
            else:
                fragment_length = 0
            lost_spacers_asstr = []
            j = fragment_length
            k = 0
            while j > 0:
                if non_zero_parent_spacers[k]:
                    loss_mask[i + k] = 0
                    lost_spacers_asstr.append(str(int(parent_spacers[i + k])))
                    j -= 1
                k += 1

            if not lost_spacers_asstr:
                self.block_loss_dict[c.name] = self.block_loss_dict.get(c.name, []) + [lost_spacers_asstr]

            i += max(k + 1, 1)
        self.block_loss_dict[c.name] = self.block_loss_dict.get(c.name, [])
        return loss_mask

    def simulation_on_tree(self, start_spacers, seed=None):
        if seed is not None:
            np.random.seed(seed)

        if self.sim_deletion_model != 'independent':
            self.block_loss_dict[self.root.name] = []

        if self.event_sim:
            if self.sim_deletion_model == 'independent':
                self.sample_events_on_tree_ilm(start_spacers)
            elif self.sim_deletion_model == 'block':
                self.sample_events_on_tree_blm(start_spacers)
            elif self.sim_deletion_model == 'fragment':
                self.sample_events_on_tree_flm(start_spacers)
            else:
                raise NotImplementedError('Simulation type not implemented!')
        else:
            if self.sim_deletion_model == 'independent':
                self.sample_spacers_on_tree_ilm(start_spacers)
            elif self.sim_deletion_model == 'block':
                self.sample_spacers_on_tree_blm(start_spacers)
            elif self.sim_deletion_model == 'fragment':
                self.sample_spacers_on_tree_flm(start_spacers)
            else:
                raise NotImplementedError('Simulation type not implemented!')
        return

    def sample_events_on_tree_ilm(self, start_spacers):
        if start_spacers == 'sample':
            start_spacers = self.sample_stationary_length_model()
        self.root.spacers = self.generate_new_spacers(self.root.name, start_spacers)
        for node in self.root.find_clades(order='preorder'):
            node.gain_loss_events = []
            parent = node.up
            self.loss_dict[node.name] = []
            self.gain_loss_events_dict[node.name] = []
            if parent is None:
                continue
            else:
                self.gain_dict[node.name] = []
            t_c = 0
            bl = node.branch_length
            spacers = copy.deepcopy(parent.spacers)
            # print(ex_spacers)

            while t_c < bl:
                # this works with setting 0 for lost spacers, might want to change this
                # ex_spacers are the indices of the actual existing spacers (!= 0) in spacers
                ex_spacers = [idx for idx, j in enumerate(spacers) if j != 0]
                insertion_time = stats.np_rate_exponential(self.gain_rate)
                # print(insertion_time)
                ls_loss_rates = [self.loss_rate] * len(ex_spacers)
                total_loss_rate = sum(ls_loss_rates)
                if ls_loss_rates:
                    loss_array_pos = np.random.choice(range(len(ex_spacers)))
                    loss_idx = ex_spacers[loss_array_pos]
                    loss_time = stats.np_rate_exponential(total_loss_rate)
                else:
                    loss_array_pos = -1
                    loss_idx = -1
                    loss_time = np.infty

                t_c += min(insertion_time, loss_time)
                if t_c > bl:
                    continue
                elif insertion_time < loss_time:
                    new_spacer = self.generate_new_spacers(node.name, 1)
                    spacers = new_spacer + spacers
                    node.gain_loss_events.append(('+', t_c, [new_spacer], len(ex_spacers), [0], spacers))
                    self.gain_loss_events_dict[node.name].append(node.gain_loss_events[-1])
                else:
                    self.loss_dict[node.name] += [spacers[loss_idx]]
                    spacers[loss_idx] = 0
                    node.gain_loss_events.append(
                        ('-', t_c, [spacers[loss_idx]], len(ex_spacers), [loss_array_pos], spacers))
                    self.gain_loss_events_dict[node.name].append(node.gain_loss_events[-1])
            node.spacers = [s for s in spacers if s != 0]

    def sample_events_on_tree_blm(self, start_spacers):
        # if seed is not None:
        #     np.random.seed(seed)
        if start_spacers == 'sample':
            start_spacers = self.sample_stationary_length_model()
        self.root.spacers = self.generate_new_spacers(self.root.name, start_spacers)

        for node in self.root.find_clades(order='preorder'):
            node.gain_loss_events = []
            parent = node.up
            self.block_loss_dict[node.name] = []
            self.loss_dict[node.name] = []
            self.gain_loss_events_dict[node.name] = []
            if parent is None:
                continue
            else:
                self.gain_dict[node.name] = []
            t_c = 0
            bl = node.branch_length
            spacers = copy.deepcopy(parent.spacers)

            while t_c < bl:
                # this works with setting 0 for lost spacers, might want to change this
                # ex_spacers are the indices of the actual existing spacers (!= 0) in spacers
                ex_spacers = [idx for idx, j in enumerate(spacers) if j != 0]
                insertion_time = stats.np_rate_exponential(self.gain_rate)
                if self.blm_rate_correction_last_spacer:
                    ls_loss_rates = [self.loss_rate * self.alpha] if len(ex_spacers) > 0 else []
                    ls_loss_rates = [self.loss_rate] * max(len(ex_spacers) - 1, 0) + ls_loss_rates
                else:
                    ls_loss_rates = [self.loss_rate] * len(ex_spacers)
                total_loss_rate = sum(ls_loss_rates)

                if ls_loss_rates:
                    # loss array pos is the actual position in the array (without 0s)
                    loss_array_pos = np.random.choice(range(len(ex_spacers)),
                                                      p=[lr / total_loss_rate for lr in ls_loss_rates])
                    # loss idx is the index in spacers (which also contains 0 of lost spacers)
                    loss_idx = ex_spacers[loss_array_pos]
                    loss_time = stats.np_rate_exponential(total_loss_rate)
                else:
                    loss_array_pos = -1
                    loss_idx = -1
                    loss_time = np.infty
                t_c += min(insertion_time, loss_time)
                if t_c > bl:
                    continue
                elif insertion_time < loss_time:
                    new_spacer = self.generate_new_spacers(node.name, 1)
                    spacers = new_spacer + spacers
                    node.gain_loss_events.append(('+', t_c, [new_spacer], len(ex_spacers), [0], spacers))
                    self.gain_loss_events_dict[node.name].append(node.gain_loss_events[-1])
                else:
                    frag_length = np.random.geometric(1 / self.alpha)
                    distance_to_first_spacer = loss_array_pos + 1
                    if self.blm_two_directional:
                        if frag_length > distance_to_first_spacer:
                            continue
                    # print('frag_length', frag_length)
                    k = 0
                    i = 0
                    lost_spacers_asstr = []
                    ls_loss_array_pos = []
                    while k < frag_length and i <= loss_idx:
                        if spacers[loss_idx - i] == 0:
                            i += 1
                            continue
                        else:
                            lost_spacers_asstr.append(str(spacers[loss_idx - i]))
                            self.loss_dict[node.name] += [spacers[loss_idx - i]]
                            spacers[loss_idx - i] = 0
                            ls_loss_array_pos.append(loss_array_pos - k)
                            i += 1
                            k += 1
                    node.gain_loss_events.append(('-', t_c, lost_spacers_asstr, len(ex_spacers), ls_loss_array_pos,
                                                  spacers))
                    self.gain_loss_events_dict[node.name].append(node.gain_loss_events[-1])
                    if not lost_spacers_asstr:
                        print(i, k, frag_length, loss_idx)
                        print(spacers)
                        raise ValueError('Empty loss set should never occur at this point of the simulation!')
                    self.block_loss_dict[node.name] = self.block_loss_dict.get(node.name, []) + [lost_spacers_asstr]
            node.spacers = [s for s in spacers if s != 0]

    def sample_events_on_tree_flm(self, start_spacers):
        if start_spacers == 'sample':
            start_spacers = self.sample_stationary_length_model()
        self.root.spacers = self.generate_new_spacers(self.root.name, start_spacers)

        for node in self.root.find_clades(order='preorder'):
            node.gain_loss_events = []
            parent = node.up
            self.block_loss_dict[node.name] = []
            self.loss_dict[node.name] = []
            self.gain_loss_events_dict[node.name] = []
            if parent is None:
                continue
            else:
                self.gain_dict[node.name] = []
            t_c = 0
            bl = node.branch_length
            spacers = copy.deepcopy(parent.spacers)

            while t_c < bl:
                ex_spacers = [idx for idx, j in enumerate(spacers) if j != 0]
                insertion_time = stats.np_rate_exponential(self.gain_rate)
                if len(ex_spacers) > 0:
                    # ls_loss_time = stats.np_rate_exponential(self.loss_rate)
                    loss_time = stats.np_rate_exponential(len(ex_spacers) * (len(ex_spacers) + 1) / 2 * self.loss_rate)
                else:
                    loss_time = np.infty

                t_c += min(insertion_time, loss_time)

                if t_c > bl:
                    continue
                elif insertion_time < loss_time:
                    new_spacer = self.generate_new_spacers(node.name, 1)
                    spacers = new_spacer + spacers

                    node.gain_loss_events.append(('+', t_c, [new_spacer], len(ex_spacers), [0], spacers))
                    self.gain_loss_events_dict[node.name].append(node.gain_loss_events[-1])
                else:
                    choice_array = list(itertools.combinations_with_replacement(range(len(ex_spacers)), 2))
                    choice_array_idx = np.random.choice(len(choice_array), 1)[0]
                    loss_array_pos = choice_array[choice_array_idx]
                    # print(len(ex_spacers), len(ex_spacers)*(len(ex_spacers) + 1) / 2, loss_array_pos)
                    # print(choice_array)
                    # loss_array_pos = np.random.choice(range(len(ex_spacers)), 2, replace=True)
                    # loss idx is the index in spacers (which also contains 0 of lost spacers)
                    loss_array_pos = sorted(loss_array_pos)
                    start_loss_idx, end_loss_idx = ex_spacers[loss_array_pos[0]], ex_spacers[loss_array_pos[1]]
                    lost_spacers_asstr = []
                    ls_loss_array_pos = []
                    i = start_loss_idx
                    j = loss_array_pos[0]
                    while i <= end_loss_idx:
                        if spacers[i] == 0:
                            i += 1
                            continue
                        else:
                            lost_spacers_asstr.append(str(spacers[i]))
                            self.loss_dict[node.name] = [spacers[i]] + self.loss_dict[node.name]
                            spacers[i] = 0
                            ls_loss_array_pos.append(j)
                            j += 1
                    lost_spacers_asstr = list(reversed(lost_spacers_asstr))
                    ls_loss_array_pos = list(reversed(ls_loss_array_pos))
                    node.gain_loss_events.append(('-', t_c, lost_spacers_asstr, len(ex_spacers), ls_loss_array_pos,
                                                  spacers))
                    self.gain_loss_events_dict[node.name].append(node.gain_loss_events[-1])
                    if not lost_spacers_asstr:
                        print(i, start_loss_idx, end_loss_idx)
                        print(spacers)
                        raise ValueError('Empty loss set should never occur at this point of the simulation!')
                    self.block_loss_dict[node.name] = self.block_loss_dict.get(node.name, []) + [lost_spacers_asstr]
            node.spacers = [s for s in spacers if s != 0]

    def sample_spacers_on_tree_ilm(self, start_spacers):
        if start_spacers == 'sample':
            start_spacers = self.sample_stationary_length_model()
        self.root.spacers = self.generate_new_spacers(self.root.name, start_spacers)

        for c in self.find_clades(order='preorder'):
            parent = c.up
            if parent is None:
                self.loss_dict[c.name] = []
                continue

            nb_gains = self.sample_gain_events(c.branch_length)
            loss_mask = self.sample_loss_events_ilm(c)
            c.spacers = np.multiply(loss_mask, parent.spacers).astype(int)
            self.loss_dict[c.name] = list(reversed([x for x in parent.spacers if x not in c.spacers]))
            c.spacers = self.generate_new_spacers(c.name, nb_gains) + list(c.spacers)
            c.spacers = [s for s in c.spacers if s != 0]

    def sample_spacers_on_tree_blm(self, start_spacers, order='preorder'):
        if start_spacers == 'sample':
            start_spacers = self.sample_stationary_length_model()
        self.root.spacers = self.generate_new_spacers(self.root.name, start_spacers)

        for c in self.find_clades(order=order):
            parent = c.up
            if parent is None:
                self.loss_dict[c.name] = []
                continue

            nb_gains = self.sample_gain_events(c.branch_length)
            loss_mask = self.sample_loss_events_blm(c)
            c.spacers = np.multiply(loss_mask, parent.spacers).astype(int)
            self.loss_dict[c.name] = list(reversed([x for x in parent.spacers if x not in c.spacers]))
            c.spacers = self.generate_new_spacers(c.name, nb_gains) + list(c.spacers)
            c.spacers = [s for s in c.spacers if s != 0]

    def duplicate_adjacent_loc_spacer(self, ex_spacers):
        idx_copy_spacer = np.random.choice(list(range(len(ex_spacers))), size=1, replace=False)
        chosen_spacer = ex_spacers[idx_copy_spacer[0]]
        new_ex_spacers = copy.deepcopy(ex_spacers)
        new_ex_spacers.insert(idx_copy_spacer[0], chosen_spacer)
        return new_ex_spacers

    def rearrange_spacer(self, ex_spacers):
        idx_copy_spacer = np.random.choice(list(range(len(ex_spacers))), size=1, replace=False)
        chosen_spacer = ex_spacers[idx_copy_spacer[0]]

        new_ex_spacers = copy.deepcopy(ex_spacers)
        del new_ex_spacers[idx_copy_spacer]
        new_loc = np.random.choice(list(range(len(new_ex_spacers))), size=1, replace=False)
        new_ex_spacers.insert(new_loc[0], chosen_spacer)
        return new_ex_spacers

    def duplicate_to_rnd_loc(self, ex_spacers):
        """
        Can place them adjacent to old spacer
        :param ex_spacers:
        :return:
        """
        idx_copy_spacer = np.random.choice(list(range(len(ex_spacers))), size=1, replace=False)
        new_loc = np.random.choice(list(range(len(ex_spacers) + 1)), size=1, replace=False)
        chosen_spacer = ex_spacers[idx_copy_spacer[0]]
        new_ex_spacers = copy.deepcopy(ex_spacers)
        new_ex_spacers.insert(new_loc[0], chosen_spacer)
        return new_ex_spacers

    def reverse_spacerblock(self, ex_spacers, block_length=2):
        """
        Throws error if block length is bigger than ex_spacers
        :param ex_spacers:
        :param block_length:
        :return:
        """
        if len(ex_spacers) == 0:
            return ex_spacers
        adjacent_reversed_blocks = list(zip(*[list(enumerate(ex_spacers))[i:] for i in range(block_length)]))
        idx_adjacent_reversed_blocks = list(range(len(adjacent_reversed_blocks)))
        idx_chosen_block = np.random.choice(idx_adjacent_reversed_blocks, size=1)
        chosen_block = adjacent_reversed_blocks[idx_chosen_block[0]]
        ls_loc = []
        ls_new_spacer = []
        for i in range(len(chosen_block)):
            ls_loc.append(chosen_block[i][0])
            ls_new_spacer.append(chosen_block[-i - 1][1])
        new_ex_spacers = copy.deepcopy(ex_spacers)
        for loc, new_spacer in zip(ls_loc, ls_new_spacer):
            new_ex_spacers[loc] = new_spacer
        return new_ex_spacers

    def sample_spacers_on_tree_flm(self, start_spacers):
        raise NotImplementedError('sample_spacers_on_tree_flm not implemented!')

    def visualize_tree(self, name='default', do_show=True, path='pictures', sort_top_order='nb_leafs',
                       determine_fsizes_by_str_len=True):
        df = pd.DataFrame(columns=['sim_spacers'])
        df_changes = pd.DataFrame(columns=['sim_gains', 'sim_losses', 'sim_block_losses'])

        top_order_sorted_by_age = []
        ls_nb_leafs = []
        ls_nodes = []
        top_order = self.get_topological_order()

        if sort_top_order == 'nb_leafs':
            for node in self.root.find_clades():
                nb_leafs = len(node.get_terminals())
                ls_nb_leafs.append(nb_leafs)
                ls_nodes.append(node)
            ls_sorted_nb_leafs, ls_sorted_nodes = (list(t) for t in zip(*sorted(zip(ls_nb_leafs, ls_nodes),
                                                                                key=lambda x: x[0])))
            ls_sorted_nb_leafs, ls_sorted_nodes = list(reversed(ls_sorted_nb_leafs)), list(reversed(ls_sorted_nodes))
            node_iterator = ls_sorted_nodes  # self.root.find_clades(order=order)
        else:
            node_iterator = self.root.find_clades(order=sort_top_order)
        for node in node_iterator:
            top_order_sorted_by_age = list(reversed(
                [s for s in self.gain_dict[node.name]])) + top_order_sorted_by_age

            sim_block_losses = self.block_loss_dict.get(node.name, [])

            df_changes.loc[node.name] = [self.gain_dict[node.name],
                                         self.loss_dict[node.name],
                                         sim_block_losses]
        top_order_sorted_by_age = [s for s in top_order if
                                   s not in top_order_sorted_by_age] + top_order_sorted_by_age
        # If multiple gains happen uniqify the topological order.
        if len(top_order_sorted_by_age) > len(top_order):
            top_order_sorted_by_age = list(
                reversed(misc.uniqify_ls_order_preserving(reversed(top_order_sorted_by_age))))
        for node in self.root.find_clades(order='level'):
            new_array = model_tools.spacer_array_to_binary(node.spacers, alignment=top_order_sorted_by_age)
            df.loc[node.name] = [new_array]

        vis.visualize_tree_sim_spacers(str(self.format(format='newick')),
                                       df_real_rec_spacers=df,
                                       df_gains_losses=df_changes,
                                       name=name,
                                       do_show=do_show,
                                       path=path,
                                       joint_del=True if self.block_loss_dict else False,
                                       top_order=top_order_sorted_by_age,
                                       determine_fsizes_by_str_len=determine_fsizes_by_str_len
                                       )

    def ascii_visualization(self, do_show=True):
        """
        Might not be bad (ASCII visualization of tree + gain/loss dicts) but not necessary
        :param do_show:
        :return:
        """
        print('Underlying tree:')
        Phylo.draw_ascii(self)
        gains_losses = {}

        # branch_labels replaces confidence values along tree -> bad formatting
        # self.gain_dict: node.name: ls_gains -> node: string(gains,losses) (required for draw function labels)
        for clade in self.root.find_clades():
            gains_losses[clade] = '+ ' + ','.join(str(e) for e in self.gain_dict[clade.name]) + '; - ' + ','.join(
                str(e) for e in self.loss_dict[clade.name])
        Phylo.draw(self, branch_labels=gains_losses, do_show=do_show)

        print('node : spacer array')
        for c in self.find_clades():
            print(c, ' : ', c.spacers)
        print('Gain-dictionary: ', self.gain_dict)
        print('Loss-dictionary: ', self.loss_dict)

    def place_dup_rear_events(self):
        something = None
        block_length = 2
        considered_nodes = self.root.find_clades()
        for node in considered_nodes:
            ex_spacers = [s for s in node.spacers if s != 0]
            if something:
                ex_spacers = self.duplicate_adjacent_loc_spacer(ex_spacers)
            elif something:
                ex_spacers = self.duplicate_to_rnd_loc(ex_spacers)
            elif something:
                ex_spacers = self.rearrange_spacer(ex_spacers)
            elif something:
                ex_spacers = self.reverse_spacerblock(ex_spacers, block_length=block_length)
            node.spacers = ex_spacers
        return
