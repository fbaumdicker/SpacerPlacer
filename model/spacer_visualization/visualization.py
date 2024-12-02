import ete3
import os
import numpy as np
import palettable
import graphviz

from model.spacer_visualization import custom_ete3


def plot_order_w_graphviz(adj_matrix, label_dict, do_show=False, save_folder=None, graph_name=None,
                          file_name=None, color_dict=None):
    graph_name = graph_name if graph_name is not None else 'PSIO'
    dot = graphviz.Digraph(graph_name, filename=file_name if file_name is not None else None,
                           node_attr={'shape': 'box', 'style': 'rounded,filled', }, )
    dot.graph_attr['rankdir'] = 'LR'
    cmap = palettable.colorbrewer.qualitative.Paired_12.hex_colors
    cmap = infinite_cmap(cmap)
    colors = [next(cmap) for _ in range(len(label_dict))] if color_dict is None \
        else [color_dict[val] for val in label_dict.values()]
    for i, label in label_dict.items():
        dot.node(str(i), label=label, fillcolor=colors[i], color='black', )

    for i in range(adj_matrix.shape[0]):
        for j in range(adj_matrix.shape[1]):
            if adj_matrix[i, j] == 1:
                dot.edge(str(i), str(j))
    dot.render(filename=file_name,
               directory=save_folder, view=do_show, format='pdf')
    dot.render(filename=file_name,
               directory=save_folder, view=do_show, format='png')


def infinite_cmap(cmap):
    while True:
        for c in cmap:
            yield c


def visualize_tree_rec_spacers(newick_tree, df_rec_spacers=None, df_gains_losses=None, name='tree', path='pictures',
                               do_show=True, rec_contra_dict=None, top_order=None, alphabet=None, opacity=1.0,
                               rec_duplications_dict=None, rec_rearrangements_dict=None, rec_double_gains_dict=None,
                               rec_indep_gains_dict=None, rec_other_dup_events_dict=None,
                               fsize=6, col_w=9, row_h=13, determine_fsizes_by_str_len=True,
                               events_fsize=4, events_col_w=7, events_row_h=7,
                               indicate_joint_del=False,
                               spacer_labels_num=True, provided_numbering=None, provided_bg_colors=None,
                               dpi=90,
                               figsize=(None, None, 'px'),
                               pool_events=True,
                               show_spacer_name_row=True,
                               pool_leaf_insertions=True):
    ete_tree_rec = ete3.Tree(newick_tree, format=3)

    ts_rec = ete3.TreeStyle()
    ts_rec.show_scale = False
    ts_rec.show_leaf_name = False
    ts_rec.margin_bottom = 10
    ts_rec.min_leaf_separation = 20
    ts_rec.show_branch_length = False

    nstyle = ete3.NodeStyle()
    nstyle["size"] = 2
    nstyle["vt_line_width"] = 1
    nstyle["hz_line_width"] = 1
    # nstyle['fgcolor'] = 'blue'
    # nstyle['partition_bgcolor'] = 'blue'

    for n in ete_tree_rec.traverse():
        n.set_style(nstyle)

    if provided_numbering is not None:
        reverse_alphabet = {val: key for key, val in alphabet.items()}
        reverse_alphabet = {key: str(provided_numbering[val]) for key, val in reverse_alphabet.items()}
        num_ls_spacers = [str(provided_numbering[c]) for c in top_order]
    else:
        reverse_alphabet = {val: key for key, val in alphabet.items()}
        num_ls_spacers = [str(alphabet[c]) for c in top_order]

    ls_spacers = [str(c) for c in top_order]

    if not spacer_labels_num or provided_numbering is not None:
        ls_gains = []
        ls_losses = []
        for idx in df_gains_losses.index:
            ls_gains.append([reverse_alphabet[s]
                             for s in df_gains_losses.loc[idx, 'rec_gains']])

            if indicate_joint_del:
                new_ls = []
                for subset in df_gains_losses.loc[idx, 'rec_losses']:
                    new_ls.append({reverse_alphabet[s] for s in subset})
                ls_losses.append(new_ls)
            else:
                ls_losses.append([reverse_alphabet[s]
                                  for s in df_gains_losses.loc[idx, 'rec_losses']])

            if rec_contra_dict is not None:
                if idx in rec_contra_dict:
                    rec_contra_dict[idx] = [reverse_alphabet[s] for s in rec_contra_dict[idx]]
            if rec_duplications_dict is not None:
                if idx in rec_duplications_dict:
                    rec_duplications_dict[idx] = [reverse_alphabet[s] for s in rec_duplications_dict[idx]]

            if rec_rearrangements_dict is not None:
                if idx in rec_rearrangements_dict:
                    rec_rearrangements_dict[idx] = [reverse_alphabet[s] for s in rec_rearrangements_dict[idx]]
            if rec_double_gains_dict is not None:
                if idx in rec_double_gains_dict:
                    rec_double_gains_dict[idx] = [reverse_alphabet[s] for s in rec_double_gains_dict[idx]]
            if rec_indep_gains_dict is not None:
                if idx in rec_indep_gains_dict:
                    rec_indep_gains_dict[idx] = [reverse_alphabet[s]
                                                 for s in rec_indep_gains_dict[idx]]
            if rec_other_dup_events_dict is not None:
                if idx in rec_other_dup_events_dict:
                    rec_other_dup_events_dict[idx] = [reverse_alphabet[s]
                                                      for s in rec_other_dup_events_dict[idx]]
        df_gains_losses['rec_gains'] = ls_gains
        df_gains_losses['rec_losses'] = ls_losses
    if determine_fsizes_by_str_len:
        if spacer_labels_num:
            col_w, events_col_w = determine_fsize(num_ls_spacers, fsize=fsize, events_fsize=events_fsize)
        else:
            col_w, events_col_w = determine_fsize(ls_spacers, fsize=fsize, events_fsize=events_fsize)

    dict_bg_colors = assign_values_to_tree(ete_tree_rec, df_rec_spacers['rec_spacers'], df_gains_losses['rec_gains'],
                                           df_gains_losses['rec_losses'], rec_contra_dict=rec_contra_dict,
                                           rec_duplications_dict=rec_duplications_dict,
                                           rec_rearrangements_dict=rec_rearrangements_dict,
                                           rec_double_gains_dict=rec_double_gains_dict,
                                           rec_indep_gains_dict=rec_indep_gains_dict,
                                           rec_other_dup_events_dict=rec_other_dup_events_dict,
                                           top_order=num_ls_spacers
                                           if spacer_labels_num or provided_numbering is not None
                                           else ls_spacers,
                                           fsize=fsize, col_w=col_w, row_h=row_h,
                                           events_fsize=events_fsize, events_col_w=events_col_w,
                                           events_row_h=events_row_h,
                                           joint_del=indicate_joint_del,
                                           provided_bg_colors=provided_bg_colors,
                                           pool_events=pool_events,
                                           pool_leaf_insertions=pool_leaf_insertions)

    fg_col_dict = {name: 'Black' for name in ls_spacers}
    bg_col_dict = {name: 'Silver' for name in ls_spacers}
    num_fg_col_dict = {name: 'Black' for name in num_ls_spacers}
    cmap = infinite_cmap(palettable.colorbrewer.qualitative.Paired_12.hex_colors)
    if provided_bg_colors is None:
        if spacer_labels_num:
            num_bg_col_dict = {name: next(cmap) for name in num_ls_spacers}
        else:
            bg_col_dict = {name: next(cmap) for name in ls_spacers}
    elif not spacer_labels_num:
        bg_col_dict = {str(key): val for key, val in provided_bg_colors.items()}
        num_bg_col_dict = {str(alphabet[key]): val for key, val in provided_bg_colors.items()}
    else:
        num_bg_col_dict = provided_bg_colors

    # Handle removing consolidated leaf insertions from the visualization
    if pool_leaf_insertions:
        all_leaf_insertions = determine_all_leaf_insertions(ete_tree_rec, df_rec_spacers['rec_spacers'],
                                                            len(ls_spacers))
        if sum(all_leaf_insertions) > 3:
            vis_ls_spacers = [s for s, a in zip(ls_spacers, all_leaf_insertions) if a != 1]
            vis_num_ls_spacers = [s for s, a in zip(num_ls_spacers, all_leaf_insertions) if a != 1]
            added_placeholder = ['+', '(' + str(sum(all_leaf_insertions)) + ')', '+']
            vis_ls_spacers = added_placeholder + vis_ls_spacers
            vis_num_ls_spacers = added_placeholder + vis_num_ls_spacers
            bg_col_dict.update({name: 'Silver' for name in added_placeholder})
            fg_col_dict.update({name: 'Black' for name in added_placeholder})
            num_bg_col_dict.update({name: 'Silver' for name in added_placeholder})
            num_fg_col_dict.update({name: 'Black' for name in added_placeholder})
        else:
            vis_ls_spacers = ls_spacers
            vis_num_ls_spacers = num_ls_spacers
    else:
        vis_ls_spacers = ls_spacers
        vis_num_ls_spacers = num_ls_spacers

    if show_spacer_name_row:
        face_spacer_nb = custom_ete3.CustomSequenceFace(vis_ls_spacers, fg_colors=fg_col_dict, bg_colors=bg_col_dict,
                                                        fsize=fsize,
                                                        col_w=col_w,
                                                        row_h=row_h)
        face_spacer_nb.inner_border.color = 'Black'
        face_spacer_nb.inner_border.width = 1
        ts_rec.aligned_header.add_face(face_spacer_nb, 1)
    if spacer_labels_num:
        num_face_spacer_nb = custom_ete3.CustomSequenceFace(vis_num_ls_spacers, fg_colors=num_fg_col_dict,
                                                            bg_colors=num_bg_col_dict,
                                                            fsize=fsize,
                                                            col_w=col_w,
                                                            row_h=row_h)
        num_face_spacer_nb.inner_border.color = 'Black'
        num_face_spacer_nb.inner_border.width = 1
        ts_rec.aligned_header.add_face(num_face_spacer_nb, 1)
    # To select random colors use next(cmap)
    legend_gain = custom_ete3.CustomSequenceFace(['2'], fg_colors={'2': 'Black'}, bg_colors={'2': 'White'},
                                                 fsize=events_fsize, col_w=events_col_w, row_h=events_row_h,
                                                 shape='ellipse',
                                                 letters_w_border={'2': 'YellowGreen'}
                                                 )
    legend_gain.margin_top = 4
    legend_gain.margin_bottom = 4
    legend_gain_text = ete3.TextFace(' Acquisitions ', fsize=events_fsize)
    legend_gain_text.margin_top = 4
    legend_gain_text.margin_bottom = 4
    ts_rec.legend.add_face(legend_gain, column=0)
    ts_rec.legend.add_face(legend_gain_text, column=1)
    legend_loss = custom_ete3.CustomSequenceFace(['3'], fg_colors={'3': 'Black'}, bg_colors={'3': 'White'},
                                                 fsize=events_fsize, col_w=events_col_w, row_h=events_row_h)
    legend_loss.inner_border.color = 'DarkRed'
    legend_loss.inner_border.width = 1
    legend_loss.margin_top = 4
    legend_loss.margin_bottom = 4
    legend_loss_text = ete3.TextFace(' Deletions ', fsize=events_fsize)
    legend_loss_text.margin_top = 4
    legend_loss_text.margin_bottom = 4
    ts_rec.legend.add_face(legend_loss, column=2)
    ts_rec.legend.add_face(legend_loss_text, column=3)
    legend_contradiction = custom_ete3.CustomSequenceFace(['5'], fg_colors={'5': 'Black'}, bg_colors={'5': 'White'},
                                                          fsize=events_fsize, col_w=events_col_w, row_h=events_row_h,
                                                          shape='hexagon',
                                                          letters_w_border={'5': 'Orange'}
                                                          )
    legend_contradiction.margin_top = 4
    legend_contradiction.margin_bottom = 4
    legend_contradiction_text = ete3.TextFace(' Contradictions ', fsize=events_fsize)
    legend_contradiction_text.margin_top = 4
    legend_contradiction_text.margin_bottom = 4
    ts_rec.legend.add_face(legend_contradiction, column=4)
    ts_rec.legend.add_face(legend_contradiction_text, column=5)
    legend_duplication = custom_ete3.CustomSequenceFace(['7'], fg_colors={'7': 'Black'}, bg_colors={'7': 'White'},
                                                        fsize=events_fsize, col_w=events_col_w, row_h=events_row_h,
                                                        shape='hexagon',
                                                        letters_w_border={'7': 'Blue'}
                                                        )
    legend_duplication.margin_top = 4
    legend_duplication.margin_bottom = 4
    legend_duplication_text = ete3.TextFace(' Duplications ', fsize=events_fsize)
    legend_duplication_text.margin_top = 4
    legend_duplication_text.margin_bottom = 4
    ts_rec.legend.add_face(legend_duplication, column=6)
    ts_rec.legend.add_face(legend_duplication_text, column=7)
    legend_rearrangement = custom_ete3.CustomSequenceFace(['11'], fg_colors={'11': 'Black'},
                                                          bg_colors={'11': 'White'},
                                                          fsize=events_fsize, col_w=events_col_w, row_h=events_row_h,
                                                          shape='hexagon',
                                                          letters_w_border={'11': 'DarkViolet'}
                                                          )
    legend_rearrangement.margin_top = 4
    legend_rearrangement.margin_bottom = 4
    legend_rearrangement_text = ete3.TextFace(' Rearrangements ', fsize=events_fsize)
    legend_rearrangement_text.margin_top = 4
    legend_rearrangement_text.margin_bottom = 4
    ts_rec.legend.add_face(legend_rearrangement, column=8)
    ts_rec.legend.add_face(legend_rearrangement_text, column=9)
    legend_double_gain = custom_ete3.CustomSequenceFace(['13'], fg_colors={'13': 'Black'}, bg_colors={'13': 'White'},
                                                        fsize=events_fsize, col_w=events_col_w, row_h=events_row_h,
                                                        shape='hexagon',
                                                        letters_w_border={'13': 'Teal'}
                                                        )
    legend_double_gain.margin_top = 4
    legend_double_gain.margin_bottom = 4
    legend_double_gain_text = ete3.TextFace(' Reacquisition ', fsize=events_fsize)
    legend_double_gain_text.margin_top = 4
    legend_double_gain_text.margin_bottom = 4
    ts_rec.legend.add_face(legend_double_gain, column=10)
    ts_rec.legend.add_face(legend_double_gain_text, column=11)

    legend_indep_gain = custom_ete3.CustomSequenceFace(['17'], fg_colors={'17': 'Black'},
                                                       bg_colors={'17': 'White'},
                                                       fsize=events_fsize, col_w=events_col_w,
                                                       row_h=events_row_h,
                                                       shape='hexagon',
                                                       letters_w_border={'17': 'SaddleBrown'}
                                                       )
    legend_indep_gain.margin_top = 4
    legend_indep_gain.margin_bottom = 4
    legend_indep_gain_text = ete3.TextFace(' Ind. acquisition ', fsize=events_fsize)
    legend_indep_gain.margin_top = 4
    legend_indep_gain.margin_bottom = 4
    ts_rec.legend.add_face(legend_indep_gain, column=12)
    ts_rec.legend.add_face(legend_indep_gain_text, column=13)

    legend_other_dup_events = custom_ete3.CustomSequenceFace(['19'], fg_colors={'19': 'Black'},
                                                             bg_colors={'19': 'White'},
                                                             fsize=events_fsize, col_w=events_col_w,
                                                             row_h=events_row_h,
                                                             shape='hexagon',
                                                             letters_w_border={'19': 'Gray'}
                                                             )
    legend_other_dup_events.margin_top = 4
    legend_other_dup_events.margin_bottom = 4
    legend_other_dup_events_text = ete3.TextFace(' other type of dup. insertion', fsize=events_fsize)
    legend_other_dup_events.margin_top = 4
    legend_other_dup_events.margin_bottom = 4
    ts_rec.legend.add_face(legend_other_dup_events, column=14)
    ts_rec.legend.add_face(legend_other_dup_events_text, column=15)

    ts_rec.legend_position = 3
    addition = '_pooled_events' if pool_leaf_insertions and pool_events else ''
    save_name = os.path.join(path, name + addition)
    ete_tree_rec.render(save_name + '.pdf', tree_style=ts_rec, dpi=dpi, w=figsize[0], h=figsize[1],
                        units=figsize[2])  # dpi works here
    ete_tree_rec.render(save_name + '.png', tree_style=ts_rec, dpi=dpi, w=figsize[0], h=figsize[1],
                        units=figsize[2])  # but not here requires fitting h, w

    if do_show:
        ete_tree_rec.show(tree_style=ts_rec, name='Reconstructed ancestors')

    return dict_bg_colors


def visualize_tree_sim_spacers(newick_tree, df_real_rec_spacers=None, df_gains_losses=None, name='tree',
                               path='pictures',
                               do_show=True, joint_del=False, top_order=None,
                               fsize=6, col_w=9, row_h=13, determine_fsizes_by_str_len=True,
                               events_fsize=4, events_col_w=7, events_row_h=7,
                               ):
    ete_tree_sim = ete3.Tree(newick_tree, format=3)

    ts_sim = ete3.TreeStyle()
    # circular tree
    # ts.mode = "c"
    ts_sim.show_scale = False
    ts_sim.show_leaf_name = False
    ts_sim.show_branch_length = False

    # ts.draw_aligned_faces_as_table = True
    # ts.rotation = 45
    # ts.branch_vertical_margin = 10
    # ts.arc_start = -180
    # ts.arc_span = 180
    # Draws nodes as small red spheres of diameter equal to 10 pixels
    # nstyle = ete3.NodeStyle()
    # nstyle["shape"] = "sphere"
    # nstyle["size"] = 4
    # nstyle["fgcolor"] = "darkred"

    # Gray dashed branch lines
    # nstyle["hz_line_type"] = 1
    # nstyle["hz_line_color"] = "#cccccc"

    # ls_spacers = list(reversed(range(1, len(df_real_rec_spacers.iloc[0, 0]) + 1)))
    top_order = [str(c) for c in top_order]

    if determine_fsizes_by_str_len:
        col_w, events_col_w = determine_fsize(top_order, fsize=fsize, events_fsize=events_fsize)
    fg_col_dict = {name: 'Black' for name in top_order}
    cmap = infinite_cmap(palettable.colorbrewer.qualitative.Paired_12.hex_colors)
    bg_col_dict = {name: next(cmap) for name in top_order}

    sim_losses = df_gains_losses['sim_block_losses'] if joint_del else df_gains_losses['sim_losses']
    assign_values_to_tree(ete_tree_sim, df_real_rec_spacers['sim_spacers'], df_gains_losses['sim_gains'],
                          sim_losses, top_order=top_order, joint_del=joint_del,
                          fsize=fsize, col_w=col_w, row_h=row_h,
                          events_fsize=events_fsize, events_col_w=events_col_w, events_row_h=events_row_h)

    face_spacer_nb = custom_ete3.CustomSequenceFace(top_order, fg_colors=fg_col_dict, bg_colors=bg_col_dict,
                                                    fsize=fsize,
                                                    col_w=col_w,
                                                    row_h=row_h)
    face_spacer_nb.inner_border.color = 'Black'
    face_spacer_nb.inner_border.width = 1

    legend_gain = custom_ete3.CustomSequenceFace(['2'], fg_colors={'2': 'Black'}, bg_colors={'2': next(cmap)},
                                                 fsize=events_fsize, col_w=events_col_w, row_h=events_row_h,
                                                 shape='ellipse',
                                                 letters_w_border={'2': 'YellowGreen'})
    ts_sim.legend.add_face(legend_gain, column=0)
    ts_sim.legend.add_face(ete3.TextFace(' Insertions ', fsize=events_fsize), column=1)
    legend_loss = custom_ete3.CustomSequenceFace(['3'], fg_colors={'3': 'Black'}, bg_colors={'3': next(cmap)},
                                                 fsize=events_fsize, col_w=events_col_w, row_h=events_row_h,
                                                 )
    legend_loss.inner_border.color = 'DarkRed'
    legend_loss.inner_border.width = 1

    ts_sim.legend_position = 3

    ts_sim.legend.add_face(legend_loss, column=2)
    ts_sim.legend.add_face(ete3.TextFace(' Deletions ', fsize=events_fsize), column=3)

    ts_sim.aligned_header.add_face(face_spacer_nb, 1)
    ts_sim.title.add_face(ete3.TextFace(name + ' simulation', fsize=fsize), 0)
    ete_tree_sim.render(os.path.join(path, name + '_sim.pdf'), tree_style=ts_sim)
    ete_tree_sim.render(os.path.join(path, name + '_sim.png'), tree_style=ts_sim)

    if do_show:
        ete_tree_sim.show(tree_style=ts_sim, name='Simulated ancestors')


def assign_values_to_tree(tree, df_spacers, gain_series, loss_series, rec_contra_dict=None, bg_cmap=None,
                          fg_cmap=None, opacity=1.0, rec_duplications_dict=None, rec_rearrangements_dict=None,
                          rec_double_gains_dict=None, rec_indep_gains_dict=None, rec_other_dup_events_dict=None,
                          joint_del=False,
                          top_order=None,
                          fsize=6, col_w=9, row_h=13,
                          events_fsize=4, events_col_w=7, events_row_h=7,
                          provided_bg_colors=None,
                          pool_events=True,
                          pool_leaf_insertions=True):
    cmap = infinite_cmap(palettable.colorbrewer.qualitative.Paired_12.hex_colors)
    default_bg_color_dict = {name: next(cmap) for name in top_order} \
        if provided_bg_colors is None else provided_bg_colors
    # default_bg_color_dict['5B'] = default_bg_color_dict['5A']
    # print(default_bg_color_dict)
    for n in tree.traverse():
        n.name = n.name.strip("''")
        if rec_double_gains_dict is None:
            set_other = set()
        else:
            set_other = set(rec_contra_dict.get(n.name, [])).union(set(rec_rearrangements_dict.get(n.name, [])),
                                                                   set(rec_duplications_dict.get(n.name, [])),
                                                                   set(rec_double_gains_dict.get(n.name, [])),
                                                                   set(rec_indep_gains_dict.get(n.name, [])),
                                                                   set(rec_other_dup_events_dict.get(n.name, [])))

        # This might be very useful for nice pictures for paper/talks! Maybe stack events into columns to lessen the
        # needed horizontal space. Horizontal space is respected by ete (for some reason).
        ls_gains = [str(c) for c in gain_series[n.name] if c not in set_other]
        ls_gains = [str(c) for c in top_order if c in ls_gains]

        fg_colors = {name: 'Black' for name in top_order}
        bg_colors = default_bg_color_dict.copy()
        if pool_events:
            if len(ls_gains) > 3:
                if n.is_leaf():
                    ls_gains_vis = ['+', '(' + str(len(ls_gains)) + ')', '+']
                    bg_colors.update({'+': 'YellowGreen', '(' + str(len(ls_gains)) + ')' : 'YellowGreen'})
                    fg_colors.update({'+': 'Black', '(' + str(len(ls_gains)) + ')' : 'Black'})
                else:
                    ls_gains_vis = [ls_gains[0], ' - ', ls_gains[-1]]
                    bg_colors.update({' - ': 'YellowGreen'})
                    fg_colors.update({' - ': 'White'})
            else:
                ls_gains_vis = ls_gains
        else:
            ls_gains_vis = ls_gains

        face_gains = custom_ete3.CustomSequenceFace(ls_gains_vis,
                                                    fg_colors=fg_colors,
                                                    bg_colors=bg_colors,
                                                    fsize=events_fsize,
                                                    col_w=events_col_w,
                                                    row_h=events_row_h,
                                                    shape='pooledHalfEllipse' if pool_events and len(ls_gains_vis) > 1
                                                    else 'ellipse',
                                                    letters_w_border={str(c): 'YellowGreen' for c in
                                                                      gain_series[n.name]})
        face_gains.margin_bottom = 2
        face_gains.opacity = opacity
        n.add_face(face_gains, 1, position="branch-top")

        # face_losses = ete3.TextFace(' ' + ', '.join([str(c) for c in loss_series[n.name]]) + ' ', fgcolor='DarkRed',
        #                             fsize=8)
        if joint_del:
            for i, j_del in enumerate(loss_series[n.name]):
                if j_del:
                    ls_del = [str(c) for c in j_del]
                    if pool_events:
                        ls_del = [str(c) for c in top_order if c in ls_del]
                    face_losses = custom_ete3.CustomSequenceFace(ls_del,
                                                                 fg_colors={name: 'Black' for name in top_order},
                                                                 bg_colors=default_bg_color_dict,
                                                                 fsize=events_fsize,
                                                                 col_w=events_col_w,
                                                                 row_h=events_row_h,
                                                                 shape='pooledRect' if pool_events else None)
                    face_losses.inner_border.color = 'DarkRed'
                    face_losses.inner_border.width = 1
                    face_losses.margin_top = 2
                    n.add_face(face_losses, i + 1, position="branch-bottom")
                    face_gains.opacity = opacity
                    face_losses.opacity = opacity

        else:
            if loss_series[n.name]:
                ls_del = [str(c) for c in loss_series[n.name]]
                if pool_events:
                    ls_del = [str(c) for c in top_order if c in ls_del]
                face_losses = custom_ete3.CustomSequenceFace(ls_del,
                                                             fg_colors={name: 'Black' for name in top_order},
                                                             bg_colors=default_bg_color_dict,
                                                             fsize=events_fsize,
                                                             col_w=events_col_w,
                                                             row_h=events_row_h,
                                                             shape='pooledRect' if pool_events else None
                                                             )
                face_losses.inner_border.color = 'DarkRed'
                face_losses.inner_border.width = 1
                face_losses.margin_top = 2
                n.add_face(face_losses, 1, position="branch-bottom")
                face_losses.opacity = opacity

        # not the optimal solution to visualize this. Maybe use sonderzeichen, bold, italic, underlined or something,
        # can only be one color
        if rec_contra_dict:
            contras = [str(c) for c in rec_contra_dict.get(n.name, [])]
            if contras:
                face_contra = custom_ete3.CustomSequenceFace(contras,
                                                             fg_colors={name: 'Black' for name in top_order},
                                                             bg_colors=default_bg_color_dict,
                                                             fsize=events_fsize,
                                                             col_w=events_col_w,
                                                             row_h=events_row_h,
                                                             shape='hexagon',
                                                             letters_w_border={c: 'Orange'
                                                                               for c in contras})
                face_contra.margin_bottom = 2
                face_contra.margin_right = 2
                n.add_face(face_contra, 0, position='branch-top')
                face_contra.opacity = opacity
        if rec_duplications_dict:
            dups = [str(c) for c in rec_duplications_dict.get(n.name, [])]
            if dups:
                face_dups = custom_ete3.CustomSequenceFace(dups,
                                                           fg_colors={name: 'Black' for name in top_order},
                                                           bg_colors=default_bg_color_dict,
                                                           fsize=events_fsize,
                                                           col_w=events_col_w,
                                                           row_h=events_row_h,
                                                           shape='hexagon',
                                                           letters_w_border={c: 'Blue'
                                                                             for c in dups}
                                                           )
                face_dups.margin_top = 1
                n.add_face(face_dups, 1, position='branch-top')
                face_dups.opacity = opacity
        if rec_rearrangements_dict:
            rearrangements = [str(c) for c in rec_rearrangements_dict.get(n.name, [])]
            if rearrangements:
                face_rearrangements = custom_ete3.CustomSequenceFace(rearrangements,
                                                                     fg_colors={name: 'Black' for name in top_order},
                                                                     bg_colors=default_bg_color_dict,
                                                                     fsize=events_fsize,
                                                                     col_w=events_col_w,
                                                                     row_h=events_row_h,
                                                                     shape='hexagon',
                                                                     letters_w_border={c: 'DarkViolet'
                                                                                       for c in rearrangements})
                face_rearrangements.margin_top = 1
                n.add_face(face_rearrangements, 1, position='branch-top')
                face_rearrangements.opacity = opacity
        if rec_double_gains_dict:
            double_gains = [str(c) for c in rec_double_gains_dict.get(n.name, [])]
            if double_gains:
                face_double_gains = custom_ete3.CustomSequenceFace(double_gains,
                                                                   fg_colors={name: 'Black' for name in top_order},
                                                                   bg_colors=default_bg_color_dict,
                                                                   fsize=events_fsize,
                                                                   col_w=events_col_w,
                                                                   row_h=events_row_h,
                                                                   shape='hexagon',
                                                                   letters_w_border={c: 'Teal'
                                                                                     for c in double_gains}
                                                                   )
                face_double_gains.margin_top = 1
                n.add_face(face_double_gains, 1, position='branch-top')
                face_double_gains.opacity = opacity

        if rec_indep_gains_dict:
            indep_gains = [str(c) for c in rec_indep_gains_dict.get(n.name, [])]
            if indep_gains:
                face_indep_gains = custom_ete3.CustomSequenceFace(indep_gains,
                                                                  fg_colors={name: 'Black' for name in
                                                                             top_order},
                                                                  bg_colors=default_bg_color_dict,
                                                                  fsize=events_fsize,
                                                                  col_w=events_col_w,
                                                                  row_h=events_row_h,
                                                                  shape='hexagon',
                                                                  letters_w_border={c: 'saddlebrown'
                                                                                    for c in
                                                                                    indep_gains}
                                                                  )
                face_indep_gains.margin_top = 1
                n.add_face(face_indep_gains, 1, position='branch-top')
                face_indep_gains.opacity = opacity
        if rec_other_dup_events_dict:
            other_dup_events = [str(c) for c in rec_other_dup_events_dict.get(n.name, [])]
            if other_dup_events:
                face_other_dup_events = custom_ete3.CustomSequenceFace(other_dup_events,
                                                                       fg_colors={name: 'Black' for name in
                                                                                  top_order},
                                                                       bg_colors=default_bg_color_dict,
                                                                       fsize=events_fsize,
                                                                       col_w=events_col_w,
                                                                       row_h=events_row_h,
                                                                       shape='hexagon',
                                                                       letters_w_border={c: 'Gray'
                                                                                         for c in
                                                                                         other_dup_events}
                                                                       )
                face_other_dup_events.margin_top = 1
                n.add_face(face_other_dup_events, 1, position='branch-top')
                face_other_dup_events.opacity = opacity

    if pool_leaf_insertions:
        all_pool_leaf_insertions = determine_all_leaf_insertions(tree, df_spacers, len(top_order))

    for n in tree.iter_leaves():
        seen_spacers = list(np.zeros(len(df_spacers.loc[n.name])))
        run_n = n
        while run_n is not None:
            seen_spacers = [1 if a == 1 or b == 1 else 0 for a, b in zip(seen_spacers, df_spacers.loc[run_n.name])]
            run_n = run_n.up

        seq = [str(t) if c == 1 else 'X' for t, c in zip(top_order, df_spacers.loc[n.name])]
        seq = [c if a == 1 else '-' for c, a in zip(seq, seen_spacers)]

        # some groundwork for using colored spacer arrays, need to think about what the spacers should be labeled
        # (int or name)
        if bg_cmap is None:
            fg_colors = {name: 'Black' for name in top_order}
            fg_colors.update({'X': 'DarkRed', '-': 'White'})
            cmap = infinite_cmap(palettable.colorbrewer.qualitative.Paired_12.hex_colors)
            bg_colors = {name: next(cmap) for name in top_order} if provided_bg_colors is None else provided_bg_colors

            bg_colors = dict(bg_colors, **{'X': 'White', '-': 'White'})
            letters_w_border = {'X': 'DarkRed'}
        else:
            bg_colors = bg_cmap
            fg_colors = fg_cmap
            letters_w_border = {'X': 'DarkRed'}

        if pool_leaf_insertions:
            if n is None:
                leaf_insertions = [0]
            else:
                leaf_insertions = [1 if a == 1 and b == 0 else 0 for a, b in zip(df_spacers.loc[n.name],
                                                                                 df_spacers.loc[n.up.name])]

            cons_seq = [a for a, b in zip(seq, leaf_insertions) if b == 1]
            len_leaf_insertions = len(cons_seq)
            if len_leaf_insertions > 3:
                insertions_symbol = '(' + str(sum(leaf_insertions)) + ')' # '-(' + str(len(cons_seq)) + ')-'
                cons_seq_vis = ['+', insertions_symbol , '+']
                bg_colors.update({insertions_symbol: 'YellowGreen', '+': 'YellowGreen'})
                fg_colors.update({insertions_symbol: 'Black', '+': 'Black'})
            else:
                if sum(all_pool_leaf_insertions) > len_leaf_insertions:
                    cons_seq_vis = ['-']*(min(sum(all_pool_leaf_insertions), 3) - len_leaf_insertions)
                    cons_seq_vis = cons_seq_vis + cons_seq
                else:
                    cons_seq_vis = cons_seq
            seq = [a for a, b in zip(seq, all_pool_leaf_insertions) if b != 1]
            seq = cons_seq_vis + seq

        face = custom_ete3.CustomSequenceFace(seq, fg_colors=fg_colors,
                                              bg_colors=bg_colors, fsize=fsize, col_w=col_w, row_h=row_h,
                                              letters_w_border=letters_w_border, letter_w_border_line_style=1)
        n.add_face(face, 1, position='aligned')
        face_name = ete3.TextFace(n.name, fsize=fsize)
        face_name.margin_left = 2
        face_name.margin_right = 2
        n.add_face(face_name, 0, position='aligned')
    return default_bg_color_dict


def determine_all_leaf_insertions(tree, df_spacers, len_top_order):
    all_pool_leaf_insertions = list(np.zeros(len_top_order))
    for n in tree.iter_leaves():
        if n is None:
            continue
        else:
            all_pool_leaf_insertions = [1 if (a == 1 and b == 0) or c == 1 else 0 for a, b, c in
                                        zip(df_spacers.loc[n.name],
                                            df_spacers.loc[n.up.name],
                                            all_pool_leaf_insertions)]
    return all_pool_leaf_insertions


def determine_fsize(str_ls_spacers, fsize=6, events_fsize=4):
    max_len = max([len(x) for x in str_ls_spacers])
    col_w = (fsize + 1) * max_len + 2
    events_col_w = (events_fsize + 1) + max_len + 2
    return col_w, events_col_w
