import collections
import copy
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from matplotlib.gridspec import GridSpec
import ast
import numpy as np

BIN_WIDTH = 0.01


# should use pickle :)
def read_csv_string_ls_to_ls(path, column):
    return pd.read_csv(path, converters={column: ast.literal_eval})


def digitize(values_to_idx_in_bins, bin_bdys):
    ls_bin_idx = []
    bin_labels = ['FS/0.0']
    bdy_prev = 0.0
    for bdy in bin_bdys[1:]:
        if bdy == 1.0:
            # bin_labels.append(f'({bdy_prev:.2f}, {bdy:.2f})')
            bin_labels.append(f'0.99')
        else:
            # bin_labels.append(f'({bdy_prev:.2f}, {bdy:.2f}]')
            bin_labels.append(f'{bdy:.2f}')
        # bin_labels.append((bdy_prev + bdy) / 2)
        bdy_prev = bdy
    for x in values_to_idx_in_bins:
        if x == 0.0:
            ls_bin_idx.append(0)
        elif x == 1.0:
            ls_bin_idx.append(len(bin_bdys))
        else:
            bdy_prev = 0.0
            for i, bdy in enumerate(bin_bdys[1:]):
                if bdy_prev < x <= bdy:
                    ls_bin_idx.append(i + 1)
                    break
                bdy_prev = bdy
    bin_labels.append('LS/1.0')
    return ls_bin_idx, bin_labels


def plot_rel_pos_freq(df, save_path, nb_cols=3, bin_width=0.04, crispr_type=None):
    array_length_bins = [0, 10, 15, 25, np.infty]
    names = ['<10', '10-15', '15-25', '25+']
    bin_width = [7, 13, 15, 31, 31]
    # array_length_bins = [0, 10, 15, 20, 30, np.infty]
    # names = ['<10', '10-15', '15-20', '20-30', '30+']
    # bin_width = [7, 13, 17, 25, 31, 31]
    df_start = df.copy()
    if crispr_type is not None:
        df_start = df_start[df_start['crispr type'] == crispr_type]
    df_start['avg array length bins'] = pd.cut(df_start['avg array length'], array_length_bins, labels=names)

    names.append('all')
    dict_df_plot = {}
    dict_nb_groups_per_bin = {'all': df_start.shape[0]}

    for j, array_length_label in enumerate(names):
        if array_length_label == 'all':
            df_plot = df_start
        else:
            df_plot = df_start[df_start['avg array length bins'] == array_length_label]

        positions = []
        for ls_ex_spacers in df_plot['nb of existent spacers'].values:
            for nb in ls_ex_spacers:
                if nb > 1:
                    positions += [i / (nb - 1) for i in range(nb)]
                else:
                    positions += [1.0]
        dict_pos_counts = collections.Counter(positions)
        all_rel_loss_pos = []
        for ls_rel_loss_pos in df_plot['relative loss positions'].values:
            all_rel_loss_pos += ls_rel_loss_pos
        dict_rel_loss_counts = collections.Counter(all_rel_loss_pos)
        # print(bin_width[j])
        # bins = [-1 / bin_width[j]] + list(np.linspace(0, 1, num=bin_width[j]))  # int(1 / bin_width)
        # bins = bins[:-1] + [0.9999999999999999] + [1 + 1 / bin_width[j]]
        bins = list(np.linspace(0, 1, num=bin_width[j]))
        # print('bins', bins)
        # p_pos_counts_bins = np.digitize(list(dict_pos_counts.keys()), bins, right=True)
        # p_rel_loss_counts_bins = np.digitize(list(dict_rel_loss_counts.keys()), bins, right=True)
        p_pos_counts_bins, bin_labels = digitize(list(dict_pos_counts.keys()), bins)
        p_rel_loss_counts_bins, _ = digitize(list(dict_rel_loss_counts.keys()), bins)
        bins = bin_labels
        # print('bin_labels', bin_labels)
        # print('val', dict_pos_counts)
        # # print('val', list(dict_pos_counts.keys()))
        # print('p_pos_counts_bins', p_pos_counts_bins)
        # print('dict_pos_counts', p_pos_counts_bins)
        binned_dict_pos_counts = {}
        binned_rel_loss_counts = {}
        for b_idx, count in zip(p_pos_counts_bins, dict_pos_counts.values()):
            binned_dict_pos_counts[bins[b_idx]] = binned_dict_pos_counts.get(bins[b_idx], 0) + count
        for b_idx, count in zip(p_rel_loss_counts_bins, dict_rel_loss_counts.values()):
            binned_rel_loss_counts[bins[b_idx]] = binned_rel_loss_counts.get(bins[b_idx], 0) + count
        # print('binned counts', binned_dict_pos_counts)
        dict_freq = {p: binned_rel_loss_counts.get(p, 0) / p_count for p, p_count in binned_dict_pos_counts.items()}
        dict_freq = [{'positions': p, 'relative loss frequencies': freq} for p, freq in dict_freq.items()]
        df_freq = pd.DataFrame(dict_freq)
        sorter = bins

        df_freq.sort_values(by="positions", key=lambda column: column.map(lambda e: sorter.index(e)), inplace=True)

        dict_df_plot[array_length_label] = df_freq
        dict_nb_groups_per_bin[array_length_label] = df_plot.shape[0]

    fig, ax = plt.subplots(int(np.ceil(len(names) / nb_cols)), nb_cols, figsize=(30, 15))
    for i, (n, df_freq) in enumerate(dict_df_plot.items()):
        row = i // nb_cols
        col = i % nb_cols

        if df_freq.empty:
            continue

        palette = ['tab:orange' if (i == 0 or i == len(df_freq['positions']) - 1) else 'tab:blue'
                   for i, _ in enumerate(df_freq['positions'])]

        if n == '15-25':
            plt.figure()
            sns.barplot(data=df_freq, x='positions', y='relative loss frequencies', palette=palette)
            plt.xlabel('relative array position')
            plt.ylabel('deletion frequency')
            plt.xticks(rotation=90)

            plt.title('Avg. array length bin ' + n + ' [%s groups]' % dict_nb_groups_per_bin[n])
            # plt.title('Deletion frequencies of groups with Avg. array length bin ' + n + '[%s groups]'
            #           % dict_nb_groups_per_bin[n])
            plt.tight_layout(pad=1.0)
            plt.savefig(save_path + '15-25.pdf')

        # ax[row, col].bar(x=df_freq['positions'], height=df_freq['relative loss frequencies'],
        #                  width=1 / (bin_width[i]),
        #                  align='center')
        # color = 'tab:blue',

        sns.barplot(data=df_freq, x='positions', y='relative loss frequencies', ax=ax[row, col], palette=palette, )

        ax[row, col].tick_params(axis='x', rotation=90)
        ax[row, col].set_xlabel('relative array position')
        ax[row, col].set_ylabel('deletion frequency')
        # sns.histplot(data=df_plot[df_plot['avg array length bins'] == n], x='relative loss positions',
        #              binwidth=BIN_WIDTH, ax=ax[row, col])
        ax[row, col].set_title('Avg. array length bin ' + n + ' [%s groups]' % dict_nb_groups_per_bin[n])

    nb_of_plots = int(np.ceil(len(names) / nb_cols)) * nb_cols

    for j in range(nb_of_plots - len(names)):
        row = - j // nb_cols - 1
        col = - j % nb_cols - 1
        ax[row, col].axis('off')

    if crispr_type is not None:
        title = 'Deletion frequencies wrt. position in array for ' + crispr_type
    else:
        title = 'Deletion frequencies wrt. position in array for all CRISPR types'

    fig.suptitle(title)
    fig.tight_layout(pad=3.0)
    fig.savefig(save_path)
    fig.show()


def plot_rel_pos(df, save_path, nb_cols=3, crispr_type=None):
    df_plot = df.explode('relative loss positions', ignore_index=True)
    if crispr_type is not None:
        df_plot = df_plot[df_plot['crispr type'] == crispr_type]
    bins = [0, 10, 15, 20, 30, np.infty]
    bin_width = [5, 13, 18, 25, 30, 30]
    names = ['<10', '10-15', '15-20', '20-30', '30+']
    df_plot['avg array length bins'] = pd.cut(df_plot['avg array length'], bins, labels=names)
    # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    names.append('all')
    fig, ax = plt.subplots(len(names) // nb_cols, nb_cols, figsize=(30, 15))
    for i, n in enumerate(names):
        if n == 'all':
            df_data = df_plot
        else:
            df_data = df_plot[df_plot['avg array length bins'] == n]
        row = i // nb_cols
        col = i % nb_cols

        sns.histplot(data=df_data, x='relative loss positions',
                     binwidth=1 / bin_width[i], ax=ax[row, col])
        ax[row, col].set_title('Avg. array length bin ' + n)
    nb_of_plots = (len(names) // nb_cols) * nb_cols
    for j in range(nb_of_plots - len(names)):
        row = - j // nb_cols - 1
        col = - j % nb_cols - 1
        ax[row, col].axis('off')
    # ax[-1, -1].axis('off')
    # sns.histplot(data=df_plot, x='relative loss positions', binwidth=BIN_WIDTH, multiple='stack',
    #              hue='avg array length bins')
    if crispr_type is not None:
        title = 'Relative loss positions for ' + crispr_type
    else:
        title = 'Relative loss positions for all CRISPR types'
    plt.suptitle(title)
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_max_loss_lengths(df, save_path, crispr_types=None, limit_x=100):
    array_count = sum([len(x) for x in df['array names']])
    df_plot = df.explode('all max loss lengths', ignore_index=True)
    nb_cols = 4
    if crispr_types is not None:
        df_plot = df_plot[df_plot['crispr type'].isin(crispr_types)]
        length_crispr_types = len(crispr_types)
    else:
        length_crispr_types = 1

    # print(df_plot)
    plt.figure()
    sns.histplot(data=df_plot, x='all max loss lengths', stat='proportion', discrete=True)
    plt.xlabel(r'deletion length')
    # limit x axis to limit_x
    plt.xlim(0, limit_x)
    plt.title(f'adjacent deletion lengths distribution [{array_count} arrays]')
    plt.tight_layout(pad=3.0)
    save_name = save_path + '_all.pdf' if crispr_types is not None else save_path
    plt.savefig(save_name)

    if crispr_types is not None:
        fig, ax = plt.subplots(length_crispr_types // nb_cols, nb_cols, figsize=(40, 15))
        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for i, n in enumerate(crispr_types):
            df_data = df_plot[df_plot['crispr type'] == n]
            if df_data.empty:
                continue
            array_count = df_data.shape[0]
            row = i // nb_cols
            col = i % nb_cols
            print(df_data['all max loss lengths'])
            sns.histplot(data=df_data, x='all max loss lengths', stat='proportion', discrete=True, ax=ax[row, col])
            ax[row, col].set_title('CRISPR type ' + n + ' [%s arrays]' % array_count)
            ax[row, col].set_xlim(0, limit_x)
            plt.xlabel(r'adjacent deletion lengths')
        plt.tight_layout(pad=3.0)
        plt.savefig(save_path)
        plt.show()


def plot_max_loss_lengths_norm(df, save_path, crispr_types=None):
    array_count = sum([len(x) for x in df['array names']])
    df_plot = df.explode('all max loss lengths (normalized)', ignore_index=True)
    # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    nb_cols = 4
    if crispr_types is not None:
        df_plot = df_plot[df_plot['crispr type'].isin(crispr_types)]

    # print(df_plot)
    plt.figure()
    sns.histplot(data=df_plot, x='all max loss lengths (normalized)', stat='proportion')
    plt.xlabel(r'deletion length (normalized)')
    plt.title(f'adjacent deletion lengths distribution (normalized) [{array_count} arrays]')
    plt.tight_layout(pad=3.0)
    save_name = save_path + '_all.pdf' if crispr_types is not None else save_path
    plt.savefig(save_name)

    if crispr_types is not None:
        fig, ax = plt.subplots(len(crispr_types) // nb_cols, nb_cols, figsize=(40, 15))
        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for i, n in enumerate(crispr_types):
            df_data = df_plot[df_plot['crispr type'] == n]
            group_count = df_data.shape[0]
            row = i // nb_cols
            col = i % nb_cols
            sns.histplot(data=df_data, x='all max loss lengths (normalized)', stat='proportion',
                         ax=ax[row, col])
            plt.xlabel(r'adjacent deletion lengths (normalized)')
            ax[row, col].set_title('CRISPR type ' + n + ' [%s groups]' % group_count)

        # sns.histplot(data=df_plot, x='all max loss lengths (normalized)', stat='proportion', binwidth=BIN_WIDTH,
        #              multiple='stack',
        #              # hue='crispr type',
        #              )
        plt.tight_layout(pad=3.0)
        plt.savefig(save_path)
        plt.show()


def plot_count_nb_unique_leafs(df, save_path, crispr_types=None):
    df_plot = df.copy()
    nb_cols = 4
    if crispr_types is not None:
        df_plot = df_plot[df_plot['crispr type'].isin(crispr_types)]
        length_crispr_types = len(crispr_types)
    else:
        length_crispr_types = 1

    # print(df_plot)
    plt.figure()
    sns.histplot(data=df_plot, x='nb leafs (after combining)', discrete=True)
    plt.xlabel(r'Number of unique CRISPR arrays')
    group_count = df_plot.shape[0]
    plt.title(f'Number of unique CRISPR arrays of all CRISPR types [{group_count} groups]')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path + '_all.pdf')

    if crispr_types is not None:
        fig, ax = plt.subplots(length_crispr_types // nb_cols, nb_cols, figsize=(40, 15))
        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for i, n in enumerate(crispr_types):
            df_data = df_plot[df_plot['crispr type'] == n]
            group_count = df_data.shape[0]
            row = i // nb_cols
            col = i % nb_cols
            sns.histplot(data=df_data, x='nb leafs (after combining)', discrete=True,
                         ax=ax[row, col])
            # sns.histplot(data=df, x='nb leafs (after combining)', discrete=True)
            plt.xlabel(r'Number of unique CRISPR arrays')
            ax[row, col].set_title('CRISPR type ' + n + ' [%s groups]' % group_count)

        plt.tight_layout(pad=3.0)
        plt.savefig(save_path)
        plt.show()


def plot_count_nb_leafs(df, save_path, crispr_types=None):
    df_plot = df.copy()
    nb_cols = 4
    if crispr_types is not None:
        df_plot = df_plot[df_plot['crispr type'].isin(crispr_types)]
        length_crispr_types = len(crispr_types)
    else:
        length_crispr_types = 1

    # print(df_plot)
    plt.figure()
    sns.histplot(data=df_plot, x='nb leafs (before combining)', discrete=True)
    plt.xlabel(r'Number of CRISPR arrays')
    group_count = df_plot.shape[0]
    plt.title(f'Number of CRISPR arrays of all CRISPR types [{group_count} groups]')
    plt.tight_layout(pad=3.0)
    save_name = save_path + '_all.pdf' if crispr_types is not None else save_path
    plt.savefig(save_name)

    if crispr_types is not None:
        fig, ax = plt.subplots(length_crispr_types // nb_cols, nb_cols, figsize=(40, 15))
        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for i, n in enumerate(crispr_types):
            df_data = df_plot[df_plot['crispr type'] == n]
            group_count = df_data.shape[0]
            row = i // nb_cols
            col = i % nb_cols
            sns.histplot(data=df_data, x='nb leafs (before combining)', discrete=True,
                         ax=ax[row, col])
            # sns.histplot(data=df, x='nb leafs (after combining)', discrete=True)
            plt.xlabel(r'Number of CRISPR arrays')
            ax[row, col].set_title('CRISPR type ' + n + ' [%s groups]' % group_count)

        plt.tight_layout(pad=3.0)
        plt.savefig(save_path)
        plt.show()


def plot_avg_array_length(df, save_path, crispr_types=None):
    # sns.histplot(data=df, x='avg array length', binwidth=1, multiple='stack', hue='crispr type')
    df_plot = df.copy()
    nb_cols = 4
    if crispr_types is not None:
        df_plot = df_plot[df_plot['crispr type'].isin(crispr_types)]
        length_crispr_types = len(crispr_types)
    else:
        length_crispr_types = 1

    # print(df_plot)
    plt.figure(figsize=(15, 5))
    sns.histplot(data=df_plot, x='avg array length', discrete=True)
    group_count = df_plot.shape[0]
    plt.xlabel(f'average array length')
    plt.title(f'average array length distribution [{group_count} groups]')
    plt.tight_layout(pad=3.0)
    save_name = save_path + '_all.pdf' if crispr_types is not None else save_path
    plt.savefig(save_name)

    if crispr_types is not None:
        fig, ax = plt.subplots(length_crispr_types // nb_cols, nb_cols, figsize=(40, 15))
        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for i, n in enumerate(crispr_types):
            df_data = df_plot[df_plot['crispr type'] == n]
            group_count = df_data.shape[0]
            row = i // nb_cols
            col = i % nb_cols
            sns.histplot(data=df_data, x='avg array length', discrete=True,
                         ax=ax[row, col])
            # sns.histplot(data=df, x='nb leafs (after combining)', discrete=True)
            plt.xlabel(r'Average array length per group')
            ax[row, col].set_title('CRISPR type ' + n + ' [%s groups]' % group_count)

        plt.tight_layout(pad=3.0)
        plt.savefig(save_path)
        plt.show()


def plot_nb_unique_spacers(df, save_path, crispr_types=None):
    # sns.histplot(data=df, x='nb gained spacers', binwidth=1)
    df_plot = df.copy()
    nb_cols = 4
    if crispr_types is not None:
        df_plot = df_plot[df_plot['crispr type'].isin(crispr_types)]
        length_crispr_types = len(crispr_types)
    else:
        length_crispr_types = 1

    # print(df_plot)
    plt.figure()
    sns.histplot(data=df_plot, x='nb unique spacers', discrete=True)
    plt.xlabel(r'Average array length per group')
    group_count = df_plot.shape[0]
    plt.title(f'Number of unique spacers per group of all CRISPR types [{group_count} groups]')
    plt.tight_layout(pad=3.0)
    save_name = save_path + '_all.pdf' if crispr_types is not None else save_path
    plt.savefig(save_name)

    if crispr_types is not None:
        fig, ax = plt.subplots(length_crispr_types // nb_cols, nb_cols, figsize=(40, 15))
        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for i, n in enumerate(crispr_types):
            df_data = df_plot[df_plot['crispr type'] == n]
            group_count = df_data.shape[0]
            row = i // nb_cols
            col = i % nb_cols
            sns.histplot(data=df_data, x='nb unique spacers', discrete=True,
                         ax=ax[row, col])
            # sns.histplot(data=df, x='nb leafs (after combining)', discrete=True)
            plt.xlabel(r'Number of unique spacers per group')
            ax[row, col].set_title('CRISPR type ' + n + ' [%s groups]' % group_count)
        plt.tight_layout(pad=3.0)
        plt.savefig(save_path)
        plt.show()


def plot_for_paper_estimated_gain_rate(df, save_path, showfliers=False, order=None, x='crispr type', rotate_90=False,
                                       figsize=(10, 5), count_threshold=None, exclude_less_strains=None,
                                       stripplot_size=2, log_transform=True):
    if count_threshold is not None:
        counts = df[x].value_counts()
        df_plot = df[df.isin(counts.index[counts >= count_threshold]).values]
    else:
        df_plot = df
    if exclude_less_strains:
        df_plot = df_plot[df_plot['nb leafs (after combining)'] >= exclude_less_strains]
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    order = df_plot[x].value_counts().index if order is None else order

    df_plot['loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    df_plot['gain_rate_1'] = df_plot['loss_rate_1_times_alpha_1'] * df_plot['avg array length']
    if log_transform:
        df_plot['gain_rate_1'] = np.log10(df_plot['gain_rate_1'])

    sns.violinplot(data=df_plot, x=x, y='gain_rate_1', hue=None, ax=ax, showfliers=showfliers,
                   order=order, cut=0, inner='quartile')
    # ax[5].set_yscale('log')
    ylims = ax.get_ylim()
    sns.stripplot(data=df_plot, x=x, y='gain_rate_1', hue=None, ax=ax, order=order, jitter=True,
                  size=stripplot_size, edgecolor='gray', linewidth=1,
                  )

    # logscale ticks
    if log_transform:
        ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax.yaxis.set_ticks(tick_range)
        ax.yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                           minor=True)
    if not showfliers:
        ax.set_ylim(ylims)
    ax.set(xlabel=r'Cas Type', ylabel=r'estimated insertion rate $\hat{\theta}_B$')
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax.set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_for_poster_estimated_parameters_vs_genus(df, save_path, showfliers=False, order=None, x='group species',
                                                  rotate_90=False,
                                                  figsize=(20, 7.5), count_threshold=None,
                                                  exclude_less_strains=None,
                                                  stripplot_size=2,
                                                  remove_extreme_outliers=False, log_transform=True):
    if count_threshold is not None:
        counts = df[x].value_counts()
        df_plot = df[df.isin(counts.index[counts >= count_threshold]).values]
    else:
        df_plot = df
    if exclude_less_strains:
        df_plot = df_plot[df_plot['nb leafs (after combining)'] >= exclude_less_strains]
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]
    order = df_plot[x].value_counts().index if order is None else order
    fig = plt.figure(figsize=figsize)
    gs = GridSpec(1, 2, figure=fig, hspace=1)
    # title = fig.add_subplot(gs[0, :])
    # fig.suptitle('Estimated BDM Parameters grouped by Genus')
    ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]

    if log_transform:
        df_plot['loss_rate_1'] = np.log10(df_plot['loss_rate_1'])
        # df_plot['alpha_1'] = np.log10(df_plot['alpha_1'])

    sns.violinplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[0], showfliers=showfliers,
                   color=sns.color_palette()[1],
                   order=order, cut=0, inner='quartile')
    # sns.boxplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
    #             order=order)
    # ax[2].set_yscale('log')
    ylims = ax[0].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[0], order=order, jitter=True, size=stripplot_size,
                  color=sns.color_palette()[1],
                  edgecolor='gray', linewidth=1)

    # logscale ticks
    if log_transform:
        ax[0].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[0].yaxis.set_ticks(tick_range)
        ax[0].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)

    if not showfliers:
        ax[0].set_ylim(ylims)
    ax[0].set(xlabel='Genus',
              ylabel=r'estimated deletion rate $\hat{\rho}_B$')  # , title=r'Block loss model: initialization rate'
    labels = [item.get_text() for item in ax[0].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[0].set_xticklabels(labels, rotation=45 if rotate_90 else 0)
    # ax[2].legend(title='Test result')
    # sns.boxplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers,
    #             order=order)

    sns.violinplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[1], showfliers=showfliers,
                   order=order, cut=0, inner='quartile', color=sns.color_palette()[1])
    # ax[3].set_yscale('log')
    ylims = ax[1].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[1], order=order, jitter=True, size=stripplot_size,
                  edgecolor='gray', linewidth=1, color=sns.color_palette()[1])
    # logscale ticks
    # if log_transform:
    #     ax[1].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    #     tick_range = np.arange(np.floor(ylims[0]), ylims[1])
    #     ax[1].yaxis.set_ticks(tick_range)
    #     ax[1].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
    #                           minor=True)
    ax[1].set(xlabel='Genus',
              ylabel=r'estimated expected block length $\hat{\alpha}$')  # title=r'Block loss model: expected block length'
    if not showfliers:
        ax[1].set_ylim(ylims)
        ax[1].set_ylim((0, 15))
        ax[1].set_yticks([x for x in list(ax[1].get_yticks()) if x != 0] + [1])
    # ax[3].legend(title='Test result')
    labels = [item.get_text() for item in ax[1].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[1].set_xticklabels(labels, rotation=45 if rotate_90 else 0)

    gs.tight_layout(fig)
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_paper_test_results_vs_type(df, save_path, order=None, x='crispr type',
                                    rotate_90=False,
                                    figsize=(20, 10), count_threshold=None,
                                    exclude_less_strains=None, ):
    if count_threshold is not None:
        counts = df[x].value_counts()
        df_plot = df[df.isin(counts.index[counts >= count_threshold]).values]
    else:
        df_plot = df
    if exclude_less_strains:
        df_plot = df_plot[df_plot['nb leafs (after combining)'] >= exclude_less_strains]
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]

    fig = plt.figure(figsize=figsize)

    gs = GridSpec(1, 1, figure=fig)
    ax = [fig.add_subplot(gs[0, 0])]
    order = df_plot[x].value_counts().index if order is None else order
    sns.countplot(data=df_plot, x=x, hue='test result', ax=ax[0],
                  order=order)
    # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    ax[0].set(title='Likelihood ratio test results',
              xlabel='')
    ax[0].legend(title='Test Result', labels=['Independent Deletion Model', 'Block Deletion Model'])
    labels = [item.get_text() for item in ax[0].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[0].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_for_poster_estimated_per_spacer_loss_rate_lrt(df, save_path, showfliers=False, order=None, x='crispr type',
                                                       rotate_90=False,
                                                       figsize=(20, 10), count_threshold=None,
                                                       exclude_less_strains=None,
                                                       stripplot_size=2,
                                                       remove_extreme_outliers=False, log_transform=True):
    if count_threshold is not None:
        counts = df[x].value_counts()
        df_plot = df[df.isin(counts.index[counts >= count_threshold]).values]
    else:
        df_plot = df
    if exclude_less_strains:
        df_plot = df_plot[df_plot['nb leafs (after combining)'] >= exclude_less_strains]
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]

    fig = plt.figure(figsize=figsize)

    gs = GridSpec(2, 2, figure=fig)
    ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    title_fig = fig.add_subplot(gs[1, :])
    title_fig.set(title='Parameter estimates grouped by Cas Type')
    title_fig.axis('off')
    order = df_plot[x].value_counts().index if order is None else order
    sns.countplot(data=df_plot, x=x, hue='test result', ax=ax[0],
                  order=order)
    # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    ax[0].set(title='Likelihood ratio test results',
              xlabel='')
    ax[0].legend(title='Test Result', labels=['Independent Deletion Model', 'Block Deletion Model'])
    labels = [item.get_text() for item in ax[0].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[0].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    df_plot['loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    if log_transform:
        df_plot['loss_rate_0'] = np.log10(df_plot['loss_rate_0'])
        df_plot['loss_rate_1_times_alpha_1'] = np.log10(df_plot['loss_rate_1_times_alpha_1'])
        df_plot['loss_rate_1'] = np.log10(df_plot['loss_rate_1'])
        # df_plot['alpha_1'] = np.log10(df_plot['alpha_1'])

    df_p = df_plot.copy()
    df_p['loss_rate_0'] = df_p['loss_rate_1_times_alpha_1']
    df_p['Deletion Model'] = [r'Block $\hat{\rho}_B \cdot \hat{\alpha}$' for _ in range(df_p.shape[0])]
    df_plot['Deletion Model'] = [r'Independent $\hat{\rho}_I$' for _ in range(df_plot.shape[0])]
    df_p = pd.concat([df_p, df_plot])
    sns.violinplot(data=df_p,
                   x=x, y='loss_rate_0', hue='Deletion Model', ax=ax[1], showfliers=showfliers,
                   order=order, inner='quartile', cut=0, hue_order=[r'Independent $\hat{\rho}_I$',
                                                                    r'Block $\hat{\rho}_B \cdot \hat{\alpha}$'])
    ylims = ax[1].get_ylim()
    sns.stripplot(data=df_p, x=x, y='loss_rate_0', hue='Deletion Model', ax=ax[1], order=order, jitter=True,
                  size=stripplot_size, dodge=True, legend=None, hue_order=[r'Independent $\hat{\rho}_I$',
                                                                           r'Block $\hat{\rho}_B \cdot \hat{\alpha}$'],
                  edgecolor='gray', linewidth=1)
    # logscale ticks
    # handles, labels = ax[0].get_legend_handles_labels()
    # l = plt.legend(handles[0:2], labels[0:2])
    if log_transform:
        ax[1].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[1].yaxis.set_ticks(tick_range)
        ax[1].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)
    if not showfliers:
        ax[1].set_ylim(ylims)
    ax[1].set(xlabel='',
              ylabel=r'estimated per spacer deletion rate $\hat{\rho}_I$, $\hat{\rho}_B \cdot \hat{\alpha}$',
              title=r'Estimated deletion rate per spacer')
    labels = [item.get_text() for item in ax[0].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[1].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    sns.violinplot(data=df_p, x=x, y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
                   color=sns.color_palette()[1],
                   order=order, cut=0, inner='quartile')
    # sns.boxplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
    #             order=order)
    # ax[2].set_yscale('log')
    ylims = ax[2].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[2], order=order, jitter=True, size=stripplot_size,
                  color=sns.color_palette()[1],
                  edgecolor='gray', linewidth=1)

    # logscale ticks
    if log_transform:
        ax[2].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[2].yaxis.set_ticks(tick_range)
        ax[2].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)

    if not showfliers:
        ax[2].set_ylim(ylims)
    ax[2].set(xlabel='Cas Type',
              ylabel=r'estimated deletion rate $\hat{\rho}_B$')  # , title=r'Block loss model: initialization rate'
    labels = [item.get_text() for item in ax[2].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[2].set_xticklabels(labels, rotation=90 if rotate_90 else 0)
    # ax[2].legend(title='Test result')
    # sns.boxplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers,
    #             order=order)

    sns.violinplot(data=df_p, x=x, y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers,
                   order=order, cut=0, inner='quartile', color=sns.color_palette()[1])
    # ax[3].set_yscale('log')
    ylims = ax[3].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[3], order=order, jitter=True, size=stripplot_size,
                  edgecolor='gray', linewidth=1, color=sns.color_palette()[1])
    # logscale ticks
    # if log_transform:
    #     ax[3].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    #     tick_range = np.arange(np.floor(ylims[0]), ylims[1])
    #     ax[3].yaxis.set_ticks(tick_range)
    #     ax[3].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
    #                           minor=True)
    ax[3].set(xlabel='Cas Type',
              ylabel=r'estimated expected block length $\hat{\alpha}$')  # title=r'Block loss model: expected block length'
    if not showfliers:
        ax[3].set_ylim(ylims)
        ax[3].set_ylim((0, 30))
        ax[3].set_yticks([x for x in list(ax[3].get_yticks()) if x != 0] + [1])
    # ax[3].legend(title='Test result')
    labels = [item.get_text() for item in ax[3].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[3].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    # plt.figtext(0.76, 0.95, 'Block Deletion Model', ha='center', va='center')

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_for_paper_estimated_per_spacer_loss_rate(df, save_path, showfliers=False, order=None, x='crispr type',
                                                  rotate_90=False,
                                                  figsize=(20, 15), count_threshold=None, exclude_less_strains=None,
                                                  stripplot_size=2,
                                                  remove_extreme_outliers=False, log_transform=True):
    if count_threshold is not None:
        counts = df[x].value_counts()
        df_plot = df[df.isin(counts.index[counts >= count_threshold]).values]
    else:
        df_plot = df
    if exclude_less_strains:
        df_plot = df_plot[df_plot['nb leafs (after combining)'] >= exclude_less_strains]
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]

    fig = plt.figure(figsize=figsize)

    sns.set_context("paper", font_scale=2.5)

    gs = GridSpec(2, 2, figure=fig)
    ax = [fig.add_subplot(gs[1, :]), fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]

    order = df_plot[x].value_counts().index if order is None else order

    df_plot['loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    if log_transform:
        df_plot['loss_rate_0'] = np.log10(df_plot['loss_rate_0'])
        df_plot['loss_rate_1_times_alpha_1'] = np.log10(df_plot['loss_rate_1_times_alpha_1'])
        df_plot['loss_rate_1'] = np.log10(df_plot['loss_rate_1'])
        # df_plot['alpha_1'] = np.log10(df_plot['alpha_1'])

    df_p = df_plot.copy()
    df_p['loss_rate_0'] = df_p['loss_rate_1_times_alpha_1']
    df_p['Deletion Model'] = [r'Block $\hat{\rho}_B \cdot \hat{\alpha}$' for _ in range(df_p.shape[0])]
    df_plot['Deletion Model'] = [r'Independent $\hat{\rho}_I$' for _ in range(df_plot.shape[0])]
    df_p = pd.concat([df_p, df_plot])
    sns.violinplot(data=df_p,
                   x=x, y='loss_rate_0', hue='Deletion Model', ax=ax[0], showfliers=showfliers,
                   order=order, inner='quartile', cut=0, hue_order=[r'Independent $\hat{\rho}_I$',
                                                                    r'Block $\hat{\rho}_B \cdot \hat{\alpha}$'])
    ylims = ax[0].get_ylim()
    sns.stripplot(data=df_p, x=x, y='loss_rate_0', hue='Deletion Model', ax=ax[0], order=order, jitter=True,
                  size=stripplot_size, dodge=True, legend=None, hue_order=[r'Independent $\hat{\rho}_I$',
                                                                           r'Block $\hat{\rho}_B \cdot \hat{\alpha}$'],
                  edgecolor='gray', linewidth=1)
    # logscale ticks
    # handles, labels = ax[0].get_legend_handles_labels()
    # l = plt.legend(handles[0:2], labels[0:2])
    if log_transform:
        ax[0].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[0].yaxis.set_ticks(tick_range)
        ax[0].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)
    if not showfliers:
        ax[0].set_ylim(ylims)
    ax[0].set(xlabel='Cas Type',
              ylabel='',
              title=r'estimated per spacer deletion rate $\hat{\rho}_I$, $\hat{\rho}_B \cdot \hat{\alpha}$',
              # title=r'Independent Deletion Model'
              )
    labels = [item.get_text() for item in ax[0].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[0].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    sns.violinplot(data=df_p, x=x, y='loss_rate_1', hue=None, ax=ax[1], showfliers=showfliers,
                   color=sns.color_palette()[1],
                   order=order, cut=0, inner='quartile')
    # sns.boxplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
    #             order=order)
    # ax[2].set_yscale('log')
    ylims = ax[1].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[1], order=order, jitter=True, size=stripplot_size,
                  color=sns.color_palette()[1],
                  edgecolor='gray', linewidth=1)

    # logscale ticks
    if log_transform:
        ax[1].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[1].yaxis.set_ticks(tick_range)
        ax[1].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)

    if not showfliers:
        ax[1].set_ylim(ylims)
    ax[1].set(xlabel='Cas type',
              ylabel='',
              title=r'estimated deletion rate $\hat{\rho}_B$')  # , title=r'Block loss model: initialization rate'
    labels = [item.get_text() for item in ax[1].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[1].set_xticklabels(labels, rotation=90 if rotate_90 else 0)
    # ax[2].legend(title='Test result')
    # sns.boxplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers,
    #             order=order)

    sns.violinplot(data=df_p, x=x, y='alpha_1', hue=None, ax=ax[2], showfliers=showfliers,
                   order=order, cut=0, inner='quartile', color=sns.color_palette()[1])
    # ax[3].set_yscale('log')
    ylims = ax[2].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[2], order=order, jitter=True, size=stripplot_size,
                  edgecolor='gray', linewidth=1, color=sns.color_palette()[1])
    # logscale ticks
    # if log_transform:
    #     ax[2].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
    #     tick_range = np.arange(np.floor(ylims[0]), ylims[1])
    #     ax[2].yaxis.set_ticks(tick_range)
    #     ax[2].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
    #                           minor=True)
    ax[2].set(xlabel='Cas type',
              ylabel='',
              title=r'estimated expected block length $\hat{\alpha}$')  # title=r'Block loss model: expected block length'
    if not showfliers:
        ax[2].set_ylim(ylims)
        ax[2].set_ylim((0, 30))
        ax[2].set_yticks([x for x in list(ax[2].get_yticks()) if x != 0] + [1])
    # ax[3].legend(title='Test result')
    labels = [item.get_text() for item in ax[2].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[2].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    # plt.figtext(0.76, 0.95, 'Block deletion model', ha='center', va='center')

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)
    sns.set_context("paper", font_scale=1)


def plot_estimated_parameters(df, save_path, showfliers=False, order=None, x='crispr type', rotate_90=False,
                              figsize=(40, 5), count_threshold=None, exclude_less_strains=None, stripplot_size=2,
                              remove_extreme_outliers=False, log_transform=True):
    if count_threshold is not None:
        counts = df[x].value_counts()
        df_plot = df[df.isin(counts.index[counts >= count_threshold]).values]
    else:
        df_plot = df
    if exclude_less_strains:
        df_plot = df_plot[df_plot['nb leafs (after combining)'] >= exclude_less_strains]
    # Check for data with no loss events, they are excluded.
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]
    fig, ax = plt.subplots(1, 6, figsize=figsize)
    order = df_plot[x].value_counts().index if order is None else order
    sns.countplot(data=df_plot, x=x, hue='test result', ax=ax[0],
                  order=order)
    # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    ax[0].set(title='Test Results',
              xlabel='Cas Type')
    ax[0].legend(title='Test Result', labels=['IDM', 'BDM'])
    labels = [item.get_text() for item in ax[0].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[0].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    df_plot['loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    df_plot['gain_rate_1'] = df_plot['loss_rate_1_times_alpha_1'] * df_plot['avg array length']
    if log_transform:
        df_plot['loss_rate_0'] = np.log10(df_plot['loss_rate_0'])
        df_plot['loss_rate_1'] = np.log10(df_plot['loss_rate_1'])
        df_plot['alpha_1'] = np.log10(df_plot['alpha_1'])

        df_plot['loss_rate_1_times_alpha_1'] = np.log10(df_plot['loss_rate_1_times_alpha_1'])
        df_plot['gain_rate_1'] = np.log10(df_plot['gain_rate_1'])

    # sns.boxplot(data=df_plot, x=x, y='loss_rate_0', hue=None, ax=ax[1], showfliers=showfliers,
    #             order=order)
    if remove_extreme_outliers:
        q = df_plot["loss_rate_0"].quantile(0.95)
        df_p = df_plot[df_plot["loss_rate_0"] < q]
    else:
        df_p = df_plot

    sns.violinplot(data=df_p,
                   x=x, y='loss_rate_0', hue=None, ax=ax[1], showfliers=showfliers,
                   order=order, inner='quartile', cut=0)
    # ax[1].set_yscale('log')
    ylims = ax[1].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='loss_rate_0', hue=None, ax=ax[1], order=order, jitter=True, size=stripplot_size,
                  edgecolor='gray', linewidth=1)
    # logscale ticks
    if log_transform:
        ax[1].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[1].yaxis.set_ticks(tick_range)
        ax[1].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)
    if not showfliers:
        ax[1].set_ylim(ylims)
    ax[1].set(xlabel='Cas Type',
              ylabel=r'estimated deletion rate $\hat{\rho}_I$', title=r'Independent Deletion Model')
    # ax[1].legend(title='Test result')
    labels = [item.get_text() for item in ax[1].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[1].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    if remove_extreme_outliers:
        q = df_plot["loss_rate_1"].quantile(0.95)
        df_p = df_plot[df_plot["loss_rate_1"] < q]
    else:
        df_p = df_plot

    sns.violinplot(data=df_p, x=x, y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
                   order=order, cut=0, inner='quartile')
    # sns.boxplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
    #             order=order)
    # ax[2].set_yscale('log')
    ylims = ax[2].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[2], order=order, jitter=True, size=stripplot_size,
                  edgecolor='gray', linewidth=1)

    # logscale ticks
    if log_transform:
        ax[2].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[2].yaxis.set_ticks(tick_range)
        ax[2].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)

    if not showfliers:
        ax[2].set_ylim(ylims)
    ax[2].set(xlabel='Cas Type',
              ylabel=r'estimated deletion rate $\hat{\rho}_B$')  # , title=r'Block loss model: initialization rate'
    labels = [item.get_text() for item in ax[2].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[2].set_xticklabels(labels, rotation=90 if rotate_90 else 0)
    # ax[2].legend(title='Test result')
    # sns.boxplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers,
    #             order=order)
    if remove_extreme_outliers:
        q = df_plot["alpha_1"].quantile(0.95)
        df_p = df_plot[df_plot["alpha_1"] < q]
    else:
        df_p = df_plot
    sns.violinplot(data=df_p, x=x, y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers,
                   order=order, cut=0, inner='quartile')
    # ax[3].set_yscale('log')
    ylims = ax[3].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='alpha_1', hue=None, ax=ax[3], order=order, jitter=True, size=stripplot_size,
                  edgecolor='gray', linewidth=1)
    # logscale ticks
    if log_transform:
        ax[3].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[3].yaxis.set_ticks(tick_range)
        ax[3].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)
    ax[3].set(xlabel='Cas Type',
              ylabel=r'estimated expected block length $\hat{\alpha}$')  # title=r'Block loss model: expected block length'
    if not showfliers:
        ax[3].set_ylim(ylims)
    # ax[3].legend(title='Test result')
    labels = [item.get_text() for item in ax[3].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[3].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    plt.figtext(0.76, 0.95, 'Block Deletion Model', ha='center', va='center')

    # sns.boxplot(data=df_plot, x=x, y='loss_rate_1_times_alpha_1',
    #             hue=None, ax=ax[4], showfliers=showfliers, order=order)
    if remove_extreme_outliers:
        q = df_plot["loss_rate_1_times_alpha_1"].quantile(0.99)
        df_p = df_plot[df_plot["loss_rate_1_times_alpha_1"] < q]
    else:
        df_p = df_plot
    sns.violinplot(data=df_p, x=x, y='loss_rate_1_times_alpha_1',
                   hue=None, ax=ax[4], showfliers=showfliers, order=order, cut=0,
                   inner='quartile')
    # ax[4].set_yscale('log')
    ylims = ax[4].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='loss_rate_1_times_alpha_1', hue=None, ax=ax[4], order=order, jitter=True,
                  size=stripplot_size, edgecolor='gray', linewidth=1)
    # logscale ticks
    if log_transform:
        ax[4].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[4].yaxis.set_ticks(tick_range)
        ax[4].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)
    if not showfliers:
        ax[4].set_ylim(ylims)
    ax[4].set(xlabel=r'Cas Type',
              ylabel=r'estimated deletion rate per spacer $\hat{\rho}_B \cdot \hat{\alpha}$')
    labels = [item.get_text() for item in ax[4].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[4].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    #  sns.boxplot(data=df_plot, x=x, y='gain_rate_1', hue=None, ax=ax[5], showfliers=showfliers,
    #              order=order)
    if remove_extreme_outliers:
        q = df_plot["gain_rate_1"].quantile(0.95)
        df_p = df_plot[df_plot["gain_rate_1"] < q]
    else:
        df_p = df_plot
    sns.violinplot(data=df_p, x=x, y='gain_rate_1', hue=None, ax=ax[5], showfliers=showfliers,
                   order=order, cut=0, inner='quartile')
    # ax[5].set_yscale('log')
    ylims = ax[5].get_ylim()
    sns.stripplot(data=df_plot, x=x, y='gain_rate_1', hue=None, ax=ax[5], order=order, jitter=True,
                  size=stripplot_size, edgecolor='gray', linewidth=1,
                  )

    # logscale ticks
    if log_transform:
        ax[5].yaxis.set_major_formatter(mticker.StrMethodFormatter("$10^{{{x:.0f}}}$"))
        tick_range = np.arange(np.floor(ylims[0]), ylims[1])
        ax[5].yaxis.set_ticks(tick_range)
        ax[5].yaxis.set_ticks([np.log10(x) for p in tick_range for x in np.linspace(10 ** p, 10 ** (p + 1), 10)],
                              minor=True)
    if not showfliers:
        ax[5].set_ylim(ylims)
    ax[5].set(xlabel=r'Cas Type', ylabel=r'estimated insertion rate $\hat{\theta}_B$')
    labels = [item.get_text() for item in ax[5].get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax[5].set_xticklabels(labels, rotation=90 if rotate_90 else 0)

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_type_counts(df, save_path, nb_arrays=False):
    # plt.figure(figsize=(20, 5))
    if nb_arrays:
        count = sum([len(x) for x in df['species by array'].tolist()])
        new_df = df.explode('species by array', ignore_index=True)
        order = new_df['crispr type'].value_counts()
    else:
        count = df.shape[0]
        order = df['crispr type'].value_counts()
    fig, ax = plt.subplots(figsize=(20, 5))
    palette = ['tab:blue' for _ in range(len(order))]
    sns.barplot(x=order.index, y=order.values, palette=palette)
    # sns.histplot(df, x='crispr type', order=order)
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels = [la.replace('Type', '') for la in labels]
    labels = [la.replace('CAS-', '') for la in labels]
    ax.set_xticklabels(labels, rotation=0)
    plt.xlabel('Cas type of array' if nb_arrays else 'Cas type of group')
    plt.ylabel('count')
    plt.title(f'Cas type distribution [{count} arrays]' if nb_arrays else f'Cas type distribution [{count} groups]')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_individual_species_counts(df, save_path):
    df_plot = df.explode('species by array', ignore_index=True)
    order = df_plot['species by array'].value_counts()
    plt.figure(figsize=(30, 5))
    palette = ['tab:blue' for _ in range(len(order))]
    sns.barplot(x=order.index, y=order.values, palette=palette)
    # sns.histplot(df_plot, x='species by array', order=order)
    plt.xticks(rotation=90)
    plt.xlabel('Species')
    plt.ylabel('Count')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_individual_genus_counts(df, save_path):
    count = sum([len(x) for x in df['genus by array'].tolist()])
    df_plot = df.explode('genus by array', ignore_index=True)
    df_plot['genus by array'] = [x.split(' ')[0] for x in df_plot['genus by array'].tolist()]
    order = df_plot['genus by array'].value_counts()
    # order = df['group species'].value_counts()
    plt.figure(figsize=(30, 10))
    # sns.histplot(df_plot, x='genus by array', order=order)
    palette = ['tab:blue' for _ in range(len(order))]
    sns.barplot(x=order.index, y=order.values, palette=palette)
    plt.xticks(rotation=90)
    plt.xlabel('genus of array')
    plt.ylabel('count')
    plt.title(f'genus distribution [{count} arrays]')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_group_species_counts(df, save_path):
    count = df.shape[0]
    plt.figure(figsize=(30, 10))
    order = df['group species'].value_counts()
    # sns.histplot(df, x='group species', order=order)
    palette = ['tab:blue' for _ in range(len(order))]
    sns.barplot(x=order.index, y=order.values, palette=palette)
    plt.xticks(rotation=90)
    plt.xlabel('genus of group')
    plt.ylabel('count')
    plt.title(f'genus distribution [{count} groups]')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_group_kingdom_counts(df, save_path):
    # print(df['kingdom'])
    # print(df['group species'])
    plt.figure(figsize=(20, 10))
    order = df['kingdom'].value_counts()
    # sns.histplot(df, x='kingdom', order=order)
    palette = ['tab:blue' for _ in range(len(order))]
    sns.barplot(x=order.index, y=order.values, palette=palette)
    plt.xlabel('Kingdom of group')
    plt.xticks(rotation=0)
    plt.ylabel('Count')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_count_orientation(df, save_path):
    df_plot = df.explode('array orientation', ignore_index=True)
    dict_keys = df_plot['array orientation'].iloc[0].keys()
    df_orient = pd.DataFrame(columns=list(dict_keys))
    orientations = df_plot['array orientation']
    for key in dict_keys:
        new_ls = []
        for o in orientations:
            new_ls.append(o[key])
        df_orient[key] = new_ls
    # extract data from dict to individual columns
    # df_plot = df_plot.join(pd.json_normalize(df['array orientation']))
    # print(df_plot)
    fig, ax = plt.subplots(1, 2, figsize=(20, 5))
    sns.countplot(data=df_orient, x='direction_prediction', ax=ax[0])
    sns.countplot(data=df_orient, x='strand_prediction', hue='strand_confidence', ax=ax[1])
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()

    # df_plot = df.explode('array orientation', ignore_index=True)
    # plt.figure(figsize=(20, 10))
    # sns.countplot(df_plot['array orientation'].astype(str), x='array orientation')
    # plt.tight_layout(pad=3.0)
    # plt.savefig(save_path)
    # plt.show()


def stack_countplot(ax):
    bottoms = {}
    for container in reversed(ax.containers):
        for bar in container:
            h = bar.get_height()
            if not np.isnan(h) and h > 0:
                x = bar.get_x()
                if x in bottoms:
                    bar.set_y(bottoms[x])
                    bottoms[x] += h
                else:
                    bottoms[x] = h
    ax.relim()
    ax.autoscale()  # recalculates the ylims due to the changed bars
    ax.yaxis.major.locator.set_params(integer=True)
    return


def plot_poster_array_orientation_agreement(df, save_path, by_group=False, filter_crispr_types=None, figsize=(10, 6)):
    if by_group:
        df_plot = copy.copy(df)
        dict_keys = df_plot['array orientation'].iloc[0][0].keys()
        # print(dict_keys)
        # print(df_plot['array orientation'].iloc[0])
        df_orient = pd.DataFrame(index=df_plot.index, columns=list(dict_keys))

        def group_orientation(ls):
            go = 0
            if all([val == 3 for val in ls]):
                go = 3
            elif 1 in set(ls) and 2 in set(ls):
                go = 3
            elif all([val in {1, 3} for val in ls]):
                go = 1
            elif all([val in {2, 3} for val in ls]):
                go = 2
            return go

        for key in dict_keys:
            new_ls = []
            for o in df_plot['array orientation'].values:
                work_ls = [val[key] for val in o]
                value = group_orientation(work_ls) if key == 'direction_prediction' else work_ls[0]
                new_ls.append(value)
            df_orient[key] = new_ls
        # print('df_orient', df_orient)
        # print(df_orient)
        df_plot = df_plot.join(df_orient)
        # print('after join', df_plot)
        # df_plot['array orientation'] = [group_orientation(ls) for ls in df['array orientation']]
    else:
        df_plot = df.explode('array orientation', ignore_index=True)
        # if isinstance(df_plot['array orientation'].iloc[0], int):
        #     dict_keys = {'direction_prediction': None, 'strand_prediction': None, 'strand_confidence': None}
        # else:
        dict_keys = df_plot['array orientation'].iloc[0].keys()
        df_orient = pd.DataFrame(index=df_plot.index, columns=list(dict_keys))
        for key in dict_keys:
            new_ls = []
            for o in df_plot['array orientation']:
                new_ls.append(o[key])
            df_orient[key] = new_ls
        # print('df_orient idx', list(df_orient.index))
        # print('df_plot idx', list(df_plot.index))
        df_plot = df_plot.join(df_orient)

    if filter_crispr_types is not None:
        # print(df_plot['crispr type'].isin(filter_crispr_types))
        # ls_ctypes = []
        # for idx in df_plot.index:
        #     ctype = df_plot.loc[idx, 'crispr type']
        #     if ctype in filter_crispr_types:
        #         ls_ctypes.append(ctype)
        #     else:
        #         ls_ctypes.append('Other')
        # df_plot['crispr type'] = ls_ctypes
        df_plot = df_plot[df_plot['crispr type'].isin(filter_crispr_types)]

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ls_we_able_to_predict = []
    ls_compared_to_other = []
    ls_other_not_confident_pred = []

    for idx in df_plot.index:
        rra = df_plot.loc[idx, 'recommend reversing array']
        sc = df_plot.loc[idx, 'strand_confidence']
        dc = df_plot.loc[idx, 'direction_prediction']
        if rra == 'nd':
            ls_we_able_to_predict.append('not confident')
        else:
            ls_we_able_to_predict.append('confident')
        # print('rra', rra)
        # print('sc', sc)
        # print('dc', dc)
        if sc != 'High' and dc == 3:
            ls_other_not_confident_pred.append('Direction & Strand not confident')
        elif sc != 'High':
            ls_other_not_confident_pred.append('Strand not confident')
        elif dc == 3:
            ls_other_not_confident_pred.append('Direction not confident')
        else:
            ls_other_not_confident_pred.append('Direction & Strand confident')

        if rra != 'nd' and sc == 'High' and dc != 3:
            ls_compared_to_other.append('All confident to predict')
        elif rra != 'nd' and sc == 'High':
            ls_compared_to_other.append('Ours + Strand')
        elif rra != 'nd' and dc != 3:
            ls_compared_to_other.append('Ours + Direction')
        elif rra != 'nd':
            ls_compared_to_other.append('Only Ours')
        else:
            ls_compared_to_other.append('Ours not confident')

    order_confidence = ['Direction & Strand not confident', 'Direction not confident', 'Strand not confident',
                        'Direction & Strand confident']
    df_plot['Prediction Confidence of Our Method'] = ls_we_able_to_predict
    df_plot['compared to Direction + Strand'] = ls_compared_to_other
    df_plot['Prediction Confidence of Direction &/or Strand'] = ls_other_not_confident_pred
    # cp = sns.countplot(data=df_plot, x='Prediction Confidence of Our Method', ax=ax[1],
    #                    hue='Prediction Confidence of Direction &/or Strand', hue_order=order_confidence, dodge=False)
    #
    # # HATCHING patches goes through like x[0], hue[0], x[1], hue[0], ....
    # hatches = ['//', '//', '//', '//', '//', '//', '', '']
    # for i, bar in enumerate(cp.patches):
    #     bar.set_hatch(hatches[i])
    #
    # stack_countplot(ax[1])
    # ax[1].legend(title='Prediction Confidence of Direction & Strand')
    # ax[1].set(title='Prediction Confidence of Direction, Strand & Our Method')
    ##############################################################################

    # filter_crispr_types = ['CAS-TypeI-B', 'CAS-TypeI-C', 'CAS-TypeI-E', 'CAS-TypeI-F',
    #                        'CAS-TypeII-A', 'CAS-TypeII-C',
    #                        'CAS-TypeIII-A', 'CAS-TypeIII-B']
    if filter_crispr_types is not None:
        # print(df_plot['crispr type'].isin(filter_crispr_types))
        # ls_ctypes = []
        # for idx in df_plot.index:
        #     ctype = df_plot.loc[idx, 'crispr type']
        #     if ctype in filter_crispr_types:
        #         ls_ctypes.append(ctype)
        #     else:
        #         ls_ctypes.append('Other')
        # df_plot['crispr type'] = ls_ctypes
        df_plot = df_plot[df_plot['crispr type'].isin(filter_crispr_types)]

    # labels = [item.get_text() for item in ax[0].legend.get_labels()]
    # labels = [la.replace('Type', '') for la in labels]
    # labels = [la.replace('CAS-', '') for la in labels]
    df_plot['crispr type'] = [s.replace('Type', '') for s in df_plot['crispr type'].values]
    df_plot['crispr type'] = [s.replace('CAS-', '') for s in df_plot['crispr type'].values]
    if filter_crispr_types is None:
        order = df_plot['crispr type'].value_counts().index
    else:
        order = [s.replace('Type', '') for s in filter_crispr_types]
        order = [s.replace('CAS-', '') for s in order]

    df_plot['direction_prediction'] = df_plot['direction_prediction'].replace({1: 'Forward',
                                                                               2: 'Reverse',
                                                                               3: 'nd'})
    our_prediction = []
    for boolean, direction in zip(df_plot['recommend reversing array'],
                                  df_plot['direction_prediction']):
        if boolean:
            op = 'Reverse' if direction == 'Forward' else 'Forward'
        else:
            op = direction
        our_prediction.append(op)

    df_plot['our_prediction'] = our_prediction

    count_all_overlap = 0
    count_direct_strand = 0
    count_direct_evo = 0
    count_strand_evo = 0
    rest = 0
    ls_states = []
    for dp, sp, op in zip(df_plot['direction_prediction'],
                          df_plot['strand_prediction'],
                          df_plot['our_prediction']):
        if dp == sp == op:
            ls_states.append('All agree')
            count_all_overlap += 1
        elif dp == sp:
            ls_states.append('Our Method disagrees')
            count_direct_strand += 1
        elif dp == op:
            ls_states.append('Strand disagrees')
            count_direct_evo += 1
        elif sp == op:
            ls_states.append('Direction disagrees')
            count_strand_evo += 1
        else:
            ls_states.append('No agreement')
            rest += 1
    df_plot['Prediction Agreement'] = ls_states
    # print('df_plot_preconfidence', df_plot)
    df_high_confidence = df_plot[df_plot['strand_confidence'] == 'High']
    df_high_confidence = df_high_confidence[df_high_confidence['recommend reversing array'] != 'nd']
    df_high_confidence = df_high_confidence[df_high_confidence['direction_prediction'] != 'nd']

    order_pa = ['All agree', 'Strand disagrees', 'Direction disagrees', 'Our Method disagrees']
    sns.countplot(data=df_high_confidence, x='Prediction Agreement', ax=ax,
                  order=order_pa, hue='crispr type', hue_order=order, dodge=True)
    stack_countplot(ax)
    ax.set(title='High Confidence Prediction Agreement btw. Direction, Strand & Our Method')

    ax.legend(title='Cas Type', fontsize=16)
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)
    return


def plot_can_predict(df, save_path, by_group=False, filter_crispr_types=None, figsize=(15, 6)):
    if by_group:
        df_plot = copy.copy(df)
        dict_keys = df_plot['array orientation'].iloc[0][0].keys()
        # print(dict_keys)
        # print(df_plot['array orientation'].iloc[0])
        df_orient = pd.DataFrame(index=df_plot.index, columns=list(dict_keys))

        def group_orientation(ls):
            go = 0
            if all([val == 3 for val in ls]):
                go = 3
            elif 1 in set(ls) and 2 in set(ls):
                go = 3
            elif all([val in {1, 3} for val in ls]):
                go = 1
            elif all([val in {2, 3} for val in ls]):
                go = 2
            return go

        for key in dict_keys:
            new_ls = []
            for o in df_plot['array orientation'].values:
                work_ls = [val[key] for val in o]
                value = group_orientation(work_ls) if key == 'direction_prediction' else work_ls[0]
                new_ls.append(value)
            df_orient[key] = new_ls
        # print('df_orient', df_orient)
        # print(df_orient)
        df_plot = df_plot.join(df_orient)
        # print('after join', df_plot)
        # df_plot['array orientation'] = [group_orientation(ls) for ls in df['array orientation']]
    else:
        df_plot = df.explode('array orientation', ignore_index=True)
        # if isinstance(df_plot['array orientation'].iloc[0], int):
        #     dict_keys = {'direction_prediction': None, 'strand_prediction': None, 'strand_confidence': None}
        # else:
        dict_keys = df_plot['array orientation'].iloc[0].keys()
        df_orient = pd.DataFrame(index=df_plot.index, columns=list(dict_keys))
        for key in dict_keys:
            new_ls = []
            for o in df_plot['array orientation']:
                new_ls.append(o[key])
            df_orient[key] = new_ls
        # print('df_orient idx', list(df_orient.index))
        # print('df_plot idx', list(df_plot.index))
        df_plot = df_plot.join(df_orient)

    if filter_crispr_types is not None:
        # print(df_plot['crispr type'].isin(filter_crispr_types))
        # ls_ctypes = []
        # for idx in df_plot.index:
        #     ctype = df_plot.loc[idx, 'crispr type']
        #     if ctype in filter_crispr_types:
        #         ls_ctypes.append(ctype)
        #     else:
        #         ls_ctypes.append('Other')
        # df_plot['crispr type'] = ls_ctypes
        df_plot = df_plot[df_plot['crispr type'].isin(filter_crispr_types)]

    fig = plt.figure(figsize=figsize)

    gs = GridSpec(1, 5, figure=fig)
    # title = fig.add_subplot(gs[0, :])
    ax = [fig.add_subplot(gs[0, :3]), fig.add_subplot(gs[0, 3:])]

    # fig, ax = plt.subplots(1, 2, figsize=(20, 10))
    ls_we_able_to_predict = []
    ls_compared_to_other = []
    ls_other_not_confident_pred = []

    for idx in df_plot.index:
        rra = df_plot.loc[idx, 'recommend reversing array']
        sc = df_plot.loc[idx, 'strand_confidence']
        dc = df_plot.loc[idx, 'direction_prediction']
        if rra == 'nd':
            ls_we_able_to_predict.append('not confident')
        else:
            ls_we_able_to_predict.append('confident')
        # print('rra', rra)
        # print('sc', sc)
        # print('dc', dc)
        if sc != 'High' and dc == 3:
            ls_other_not_confident_pred.append('Direction, Strand not confid.')
        elif sc != 'High':
            ls_other_not_confident_pred.append('Strand not confid.')
        elif dc == 3:
            ls_other_not_confident_pred.append('Direction not confid.')
        else:
            ls_other_not_confident_pred.append('Direction, Strand confid.')

        if rra != 'nd' and sc == 'High' and dc != 3:
            ls_compared_to_other.append('All confident to predict')
        elif rra != 'nd' and sc == 'High':
            ls_compared_to_other.append('Ours + Strand')
        elif rra != 'nd' and dc != 3:
            ls_compared_to_other.append('Ours + Direction')
        elif rra != 'nd':
            ls_compared_to_other.append('Only Ours')
        else:
            ls_compared_to_other.append('Ours not confident')

    order_confidence = ['Direction, Strand not confid.', 'Direction not confid.', 'Strand not confid.',
                        'Direction, Strand confid.']
    df_plot['Prediction Confidence of Our Method'] = ls_we_able_to_predict
    df_plot['compared to Direction + Strand'] = ls_compared_to_other
    df_plot['Prediction Confidence of Direction &/or Strand'] = ls_other_not_confident_pred
    cp = sns.countplot(data=df_plot, x='Prediction Confidence of Our Method', ax=ax[1],
                       hue='Prediction Confidence of Direction &/or Strand', hue_order=order_confidence, dodge=False)

    # HATCHING patches goes through like x[0], hue[0], x[1], hue[0], ....
    hatches = ['//', '//', '//', '//', '//', '//', '', '']
    for i, bar in enumerate(cp.patches):
        bar.set_hatch(hatches[i])

    stack_countplot(ax[1])
    ax[1].legend(title='Prediction Confidence')
    ax[1].set(title='Prediction Confidence of Direction, Strand & Our Method')
    ##############################################################################

    # filter_crispr_types = ['CAS-TypeI-B', 'CAS-TypeI-C', 'CAS-TypeI-E', 'CAS-TypeI-F',
    #                        'CAS-TypeII-A', 'CAS-TypeII-C',
    #                        'CAS-TypeIII-A', 'CAS-TypeIII-B']
    if filter_crispr_types is not None:
        # print(df_plot['crispr type'].isin(filter_crispr_types))
        # ls_ctypes = []
        # for idx in df_plot.index:
        #     ctype = df_plot.loc[idx, 'crispr type']
        #     if ctype in filter_crispr_types:
        #         ls_ctypes.append(ctype)
        #     else:
        #         ls_ctypes.append('Other')
        # df_plot['crispr type'] = ls_ctypes
        df_plot = df_plot[df_plot['crispr type'].isin(filter_crispr_types)]

    # labels = [item.get_text() for item in ax[0].legend.get_labels()]
    # labels = [la.replace('Type', '') for la in labels]
    # labels = [la.replace('CAS-', '') for la in labels]
    df_plot['crispr type'] = [s.replace('Type', '') for s in df_plot['crispr type'].values]
    df_plot['crispr type'] = [s.replace('CAS-', '') for s in df_plot['crispr type'].values]
    if filter_crispr_types is None:
        order = df_plot['crispr type'].value_counts().index
    else:
        order = [s.replace('Type', '') for s in filter_crispr_types]
        order = [s.replace('CAS-', '') for s in order]

    df_plot['direction_prediction'] = df_plot['direction_prediction'].replace({1: 'Forward',
                                                                               2: 'Reverse',
                                                                               3: 'nd'})
    our_prediction = []
    for boolean, direction in zip(df_plot['recommend reversing array'],
                                  df_plot['direction_prediction']):
        if boolean:
            op = 'Reverse' if direction == 'Forward' else 'Forward'
        else:
            op = direction
        our_prediction.append(op)

    df_plot['our_prediction'] = our_prediction

    count_all_overlap = 0
    count_direct_strand = 0
    count_direct_evo = 0
    count_strand_evo = 0
    rest = 0
    ls_states = []
    for dp, sp, op in zip(df_plot['direction_prediction'],
                          df_plot['strand_prediction'],
                          df_plot['our_prediction']):
        if dp == sp == op:
            ls_states.append('All agree')
            count_all_overlap += 1
        elif dp == sp:
            ls_states.append('Our Method disagrees')
            count_direct_strand += 1
        elif dp == op:
            ls_states.append('Strand disagrees')
            count_direct_evo += 1
        elif sp == op:
            ls_states.append('Direction disagrees')
            count_strand_evo += 1
        else:
            ls_states.append('No agreement')
            rest += 1
    df_plot['Prediction Agreement'] = ls_states
    # print('df_plot_preconfidence', df_plot)
    df_high_confidence = df_plot[df_plot['strand_confidence'] == 'High']
    df_high_confidence = df_high_confidence[df_high_confidence['recommend reversing array'] != 'nd']
    df_high_confidence = df_high_confidence[df_high_confidence['direction_prediction'] != 'nd']

    order_pa = ['All agree', 'Strand disagrees', 'Direction disagrees', 'Our Method disagrees']
    sns.countplot(data=df_high_confidence, x='Prediction Agreement', ax=ax[0],
                  order=order_pa, hue='crispr type', hue_order=order, dodge=True)
    stack_countplot(ax[0])
    ax[0].set(title='High Confidence Prediction Agreement btw. Direction, Strand & Our Method')

    ax[0].legend(title='Cas Type', fontsize='x-large')
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)
    return


def plot_orientation_stats(df, save_path, filter_3=False, by_group=False, filter_more_than_3_unique=False,
                           percentage=False, filter_crispr_types=None):
    # if large_crispr_types is None:
    #     large_crispr_types = []
    if by_group:
        df_plot = copy.copy(df)
        dict_keys = df_plot['array orientation'].iloc[0][0].keys()
        # print(dict_keys)
        # print(df_plot['array orientation'].iloc[0])
        df_orient = pd.DataFrame(index=df_plot.index, columns=list(dict_keys))

        def group_orientation(ls):
            go = 0
            if all([val == 3 for val in ls]):
                go = 3
            elif 1 in set(ls) and 2 in set(ls):
                go = 3
            elif all([val in {1, 3} for val in ls]):
                go = 1
            elif all([val in {2, 3} for val in ls]):
                go = 2
            return go

        for key in dict_keys:
            new_ls = []
            for o in df_plot['array orientation'].values:
                work_ls = [val[key] for val in o]
                value = group_orientation(work_ls) if key == 'direction_prediction' else work_ls[0]
                new_ls.append(value)
            df_orient[key] = new_ls
        # print('df_orient', df_orient)
        # print(df_orient)
        df_plot = df_plot.join(df_orient)
        # print('after join', df_plot)
        # df_plot['array orientation'] = [group_orientation(ls) for ls in df['array orientation']]
    else:
        df_plot = df.explode('array orientation', ignore_index=True)
        dict_keys = df_plot['array orientation'].iloc[0].keys()
        df_orient = pd.DataFrame(index=df_plot.index, columns=list(dict_keys))
        for key in dict_keys:
            new_ls = []
            for o in df_plot['array orientation']:
                new_ls.append(o[key])
            df_orient[key] = new_ls
        # print('df_orient idx', list(df_orient.index))
        # print('df_plot idx', list(df_plot.index))
        df_plot = df_plot.join(df_orient)

    if filter_3:
        df_plot = df_plot[df_plot['direction_prediction'] != 3]
    if filter_more_than_3_unique:
        df_plot = df_plot[df_plot['nb leafs (after combining)'] > 3]
    fig, ax = plt.subplots(2, 4, figsize=(30, 10))

    add_str = ' '
    if filter_crispr_types is not None:
        # print(df_plot['crispr type'].isin(filter_crispr_types))
        df_plot = df_plot[df_plot['crispr type'].isin(filter_crispr_types)]
        # print(df_plot.shape[0])
        add_str += ' '.join(filter_crispr_types) if len(filter_crispr_types) < 5 else 'many types'
    if by_group:
        if filter_more_than_3_unique:
            title = 'Orientation stats of groups with more than 3 unique arrays' + add_str
        else:
            title = 'Orientation stats of groups' + add_str
        plt.suptitle(title)
    else:
        plt.suptitle('Orientation stats of arrays' + add_str)
    order_rr = [True, False, 'nd']
    order_dp = [1, 2, 3]
    order_strand = ['Forward', 'Reverse']
    order_opostr = ['Forward', 'Reverse', 'nd']
    order_opodi = [1, 2, 'nd']
    order_pa = ['All agree', 'Direction + Strand', 'Direction + Our', 'Strand + Our']
    order_conf = ['High', 'Medium', 'Low']
    # print('recommend', df_plot.shape[0])
    # print('recommend type 2', df_plot[df_plot['crispr type'] == 'CAS-TypeII-A'].shape[0])
    sns.countplot(data=df_plot, x='recommend reversing array', ax=ax[0, 0], order=order_rr,
                  hue='crispr type', hue_order=filter_crispr_types, dodge=False)
    stack_countplot(ax[0, 0])
    # df_plot_s = df.copy()
    # df_plot_s['recommend reversing array'] = df_plot_s['recommend reversing array'].replace({'nd': 2})
    # sns.histplot(data=df_plot_s, x='recommend reversing array',
    #             ax=ax[0, 0],
    #             hue='crispr type', hue_order=filter_crispr_types, multiple='stack', discrete=True)
    # ax[0, 0].set_xticklabels(order_rr)
    ax[0, 0].set(title='Recommended reversion' + add_str)
    # ax[0].legend(title='Recommended reversion', labels=['IDM', 'BDM'])
    # print('direction', df_plot.shape[0])
    # print('direction type 2', df_plot[df_plot['crispr type'] == 'CAS-TypeII-A'].shape[0])
    sns.countplot(data=df_plot, x='direction_prediction', ax=ax[0, 1], order=order_dp, hue='crispr type',
                  hue_order=filter_crispr_types, dodge=False)
    stack_countplot(ax[0, 1])
    ax[0, 1].set(title='CRISPRDirection orientation' + add_str)
    mask = df_plot['recommend reversing array'].values

    def pred_f(x, m):
        if m == 'nd':
            return 'nd'
            # if x == 1:
            #     return 'nd (1)'
            # if x == 2:
            #     return 'nd (2)'
            # if x == 3:
            #     return 'nd (3)'
        elif m:
            if x == 1:
                return 2
            if x == 2:
                return 1
            if x == 3:
                return 2
                # return '2 (3)'
        else:
            if x == 3:
                return 1
                # return '1 (3)'
            else:
                return x

    df_plot['our predicted orientation'] = [pred_f(o, m)
                                            for o, m in zip(df_plot['direction_prediction'].values, mask)]
    sns.countplot(data=df_plot, x='direction_prediction', ax=ax[1, 1], hue='our predicted orientation',
                  order=order_dp, hue_order=order_opodi)
    ax[1, 1].set(title='our predicted orientation vs direction' + add_str)

    sns.countplot(data=df_plot, x='strand_prediction', ax=ax[0, 2], order=order_strand, hue='strand_confidence',
                  hue_order=order_conf)
    ax[0, 2].set(title='CRISPRStrand orientation' + add_str)
    mask = df_plot['recommend reversing array'].values

    def pred_f_strand(x, m):
        if m == 'nd':
            return 'nd'
            # if x == 1:
            #     return 'nd (1)'
            # if x == 2:
            #     return 'nd (2)'
            # if x == 3:
            #     return 'nd (3)'
        elif m:
            if x == 'Forward':
                return 'Reverse'
            if x == 'Reverse':
                return 'Forward'
        else:
            return x

    df_plot['our predicted orientation'] = [pred_f_strand(o, m)
                                            for o, m in zip(df_plot['strand_prediction'].values, mask)]

    df_high_confidence = df_plot[df_plot['strand_confidence'] == 'High']
    sns.countplot(data=df_high_confidence, x='strand_prediction', ax=ax[1, 2],
                  hue='our predicted orientation', order=order_strand, hue_order=order_opostr)
    ax[1, 2].set(title='our prediction orientation vs strand (high confidence)' + add_str)

    sns.countplot(data=df_plot, x='strand_prediction', ax=ax[0, 3], hue='our predicted orientation',
                  order=order_strand, hue_order=order_opostr)
    ax[0, 3].set(title='our predicted orientation vs strand (all)' + add_str)

    df_low_confidence = df_plot[df_plot['strand_confidence'] != 'High']
    sns.countplot(data=df_low_confidence, x='strand_prediction', ax=ax[1, 3], hue='our predicted orientation',
                  order=order_strand, hue_order=order_opostr)
    ax[1, 3].set(title='our predicted orientation vs strand (low/middle confidence)' + add_str)

    df_plot['direction_prediction'] = df_plot['direction_prediction'].replace({1: 'Forward',
                                                                               2: 'Reverse',
                                                                               3: 'nd'})
    our_prediction = []
    for boolean, direction in zip(df_plot['recommend reversing array'],
                                  df_plot['direction_prediction']):
        if boolean:
            op = 'Reverse' if direction == 'Forward' else 'Forward'
        else:
            op = direction
        our_prediction.append(op)
    df_plot['our_prediction'] = our_prediction

    count_all_overlap = 0
    count_direct_strand = 0
    count_direct_evo = 0
    count_strand_evo = 0
    rest = 0
    ls_states = []
    for dp, sp, op in zip(df_plot['direction_prediction'],
                          df_plot['strand_prediction'],
                          df_plot['our_prediction']):
        if dp == sp == op:
            ls_states.append('All agree')
            count_all_overlap += 1
        elif dp == sp:
            ls_states.append('Direction + Strand')
            count_direct_strand += 1
        elif dp == op:
            ls_states.append('Direction + Our')
            count_direct_evo += 1
        elif sp == op:
            ls_states.append('Strand + Our')
            count_strand_evo += 1
        else:
            ls_states.append('No agreement')
            rest += 1
    df_plot['prediction agreement'] = ls_states
    # print('df_plot_preconfidence', df_plot)
    df_high_confidence = df_plot[df_plot['strand_confidence'] == 'High']
    df_high_confidence = df_high_confidence[df_high_confidence['recommend reversing array'] != 'nd']
    df_high_confidence = df_high_confidence[df_high_confidence['direction_prediction'] != 'nd']
    # print(df_high_confidence['prediction agreement'].value_counts())
    sns.countplot(data=df_high_confidence, x='prediction agreement', ax=ax[1, 0],
                  order=order_pa, hue='crispr type', hue_order=filter_crispr_types, dodge=False)
    stack_countplot(ax[1, 0])
    ax[1, 0].set(title='prediction agreement (high confidence)' + add_str)

    # df_low_confidence = df_plot[df_plot['strand_confidence'] != 'High']
    # df_low_confidence = df_low_confidence[df_low_confidence['recommend reversing array'] == 'nd']
    # df_low_confidence = df_low_confidence[df_low_confidence['direction_prediction'] == 'nd']
    #
    # sns.countplot(data=df_low_confidence, x='Prediction agreement', ax=ax[2, 0])
    # ax[2, 0].set(title='Prediction agreement (low confidence)')

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_orientation_likelihood_distribution(df, save_path, log=True, filter_3=False, filter_crispr_types=None):
    # lh_0_key = ''
    # reversed_lh_0_key = ''
    add_str = ' '
    if filter_crispr_types is not None:
        df_plot = df[df['crispr type'].isin(filter_crispr_types)]
        add_str += ' '.join(filter_crispr_types) if len(filter_crispr_types) < 5 else 'many types'
    else:
        df_plot = df.copy()

    def f(col):
        mask = []
        for i in col.values:
            if 3 not in i:
                mask.append(True)
            else:
                mask.append(False)
        return mask

    if filter_3:
        df_plot = df_plot[f(df_plot['array orientation'])]
    if log:
        lh_key = 'ln_lh_1'
        reversed_lh_key = 'reversed_ln_lh_1'
        lh_ratio = 'ln_lh_1 - reversed_ln_lh_1'
    else:
        lh_key = 'lh_1'
        reversed_lh_key = 'reversed_lh_1'
        lh_ratio = 'lh_1 / reversed_lh_1'
    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    sns.histplot(data=df_plot, x=lh_key, ax=ax[0], stat='probability')
    sns.histplot(data=df_plot, x=reversed_lh_key, ax=ax[1], stat='probability')
    sns.histplot(data=df_plot, x=lh_ratio, hue='recommend reversing array', multiple='stack', ax=ax[2],
                 stat='probability',
                 binwidth=1)
    ax[2].set_xlim(-20, 20)
    if filter_3:
        fig.suptitle('log likelihood distributions (orientation) wo 3s'
                     if log else 'likelihood distributions (orientation) wo 3s')
    else:
        fig.suptitle('log likelihood distributions (orientation)' + add_str
                     if log else 'likelihood distributions (orientation)' + add_str)
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_fes_les_loss_freq(df, save_path, max_considered_distance=10, use_global=False, crispr_type='all types'):
    if crispr_type != 'all types':
        df_plot = df[df['crispr type'] == crispr_type]
    else:
        df_plot = df.copy()
    add_str = ' ' + crispr_type

    fes_key = 'global fes +m presence' if use_global else 'fes +m presence'
    les_key = 'global les -m presence' if use_global else 'les -m presence'
    ls_key = 'global ls -m presence' if use_global else 'ls -m presence'

    ls_fes_presences = [0] * (max_considered_distance + 1)
    ls_fes_totals = [0] * (max_considered_distance + 1)
    ls_fes_labels = ['+' + str(m) for m in range(max_considered_distance + 1)]
    ls_les_presences = [0] * (max_considered_distance + 1)
    ls_les_totals = [0] * (max_considered_distance + 1)
    ls_les_labels = ['-' + str(m) for m in range(max_considered_distance + 1)]
    ls_ls_presences = [0] * (max_considered_distance + 1)
    ls_ls_totals = [0] * (max_considered_distance + 1)
    ls_ls_labels = ['-' + str(m) for m in range(max_considered_distance + 1)]

    for dict_fes_presence in df_plot[fes_key].values:
        for m, val in dict_fes_presence.items():
            if m <= max_considered_distance:
                ls_fes_presences[m] += val[0]
                ls_fes_totals[m] += val[1]

    for dict_les_presence in df_plot[les_key].values:
        for m, val in dict_les_presence.items():
            if m <= max_considered_distance:
                ls_les_presences[m] += val[0]
                ls_les_totals[m] += val[1]

    for dict_ls_presence in df_plot[ls_key].values:
        for m, val in dict_ls_presence.items():
            if m <= max_considered_distance:
                ls_ls_presences[m] += val[0]
                ls_ls_totals[m] += val[1]
    ls_fes_freq = []
    ls_les_freq = []
    ls_ls_freq = []
    for p, t in zip(ls_fes_presences, ls_fes_totals):
        if t > 0:
            ls_fes_freq.append(1 - p / t)
        else:
            ls_fes_freq.append(0)
    for p, t in zip(ls_les_presences, ls_les_totals):
        if t > 0:
            ls_les_freq.append(1 - p / t)
        else:
            ls_les_freq.append(0)
    for p, t in zip(ls_ls_presences, ls_ls_totals):
        if t > 0:
            ls_ls_freq.append(1 - p / t)
        else:
            ls_ls_freq.append(0)

    # avg_deletion_rate_1 = df['loss_rate_1'].mean() * df['alpha_1'].mean()
    # avg_deletion_rate_0 = df['loss_rate_0'].mean()

    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    fig.suptitle('Global FES/LES/LS deletion frequencies' + add_str
                 if use_global else 'FES/LES/LS deletion frequencies' + add_str)
    ax[0].bar(x=ls_fes_labels, height=ls_fes_freq, width=0.8, align='center')
    ax[0].set_xlabel('Positions to the right of FES')
    ax[0].set_ylabel('Deletion frequencies')
    ax[0].set_title('FES+ deletion frequencies' + add_str)

    ax[1].bar(x=ls_les_labels[::-1], height=ls_les_freq[::-1], width=0.8, align='center')
    ax[1].set_xlabel('Positions to the left of LES')
    ax[1].set_ylabel('Deletion frequencies')
    ax[1].set_title('LES- deletion frequencies' + add_str)

    if not use_global:
        ax[2].bar(x=ls_ls_labels[::-1], height=ls_ls_freq[::-1], width=0.8, align='center')
        ax[2].set_xlabel('Positions to the left of LS')
        ax[2].set_ylabel('Deletion frequencies')
        ax[2].set_title('LS- deletion frequencies' + add_str)

        ymin = min(ax[0].get_ylim()[0], ax[1].get_ylim()[0], ax[2].get_ylim()[0])
        ymax = max(ax[0].get_ylim()[1], ax[1].get_ylim()[1], ax[2].get_ylim()[1])
    else:
        ax[2].axis('off')
        ymin = min(ax[0].get_ylim()[0], ax[1].get_ylim()[0])
        ymax = max(ax[0].get_ylim()[1], ax[1].get_ylim()[1])

    ax[0].set_ylim([ymin, ymax])
    ax[1].set_ylim([ymin, ymax])
    ax[2].set_ylim([ymin, ymax])

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_paper_fs_ls_loss_freq(df, save_path, max_considered_distance=10, binning=False, interlaced=False):
    df_plot = df.copy()
    sns.set_context("paper", font_scale=2)
    if binning:
        bins = [0, 15, 25, np.infty]
        names = ['<15', '15 - 25', '>25']
        df_plot['avg array length bins'] = pd.cut(df_plot['avg array length'], bins, labels=names)
        df_plot = df_plot[df_plot['avg array length bins'] == '15 - 25']

    fs_key = 'fs +m presence'
    ls_key = 'ls -m presence'

    ls_fs_presences = [0] * (max_considered_distance + 1)
    ls_fs_totals = [0] * (max_considered_distance + 1)
    ls_fs_labels = ['FS'] + ['+' + str(m) for m in range(1, max_considered_distance + 1)]
    ls_ls_presences = [0] * (max_considered_distance + 1)
    ls_ls_totals = [0] * (max_considered_distance + 1)
    ls_ls_labels = ['LS'] + ['-' + str(m) for m in range(1, max_considered_distance + 1)]

    for dict_fs_presence in df_plot[fs_key].values:
        for m, val in dict_fs_presence.items():
            if m <= max_considered_distance:
                ls_fs_presences[m] += val[0]
                ls_fs_totals[m] += val[1]

    for dict_ls_presence in df_plot[ls_key].values:
        for m, val in dict_ls_presence.items():
            if m <= max_considered_distance:
                ls_ls_presences[m] += val[0]
                ls_ls_totals[m] += val[1]
    ls_fs_freq = []
    ls_ls_freq = []
    for p, t in zip(ls_fs_presences, ls_fs_totals):
        if t > 0:
            ls_fs_freq.append(1 - p / t)
        else:
            ls_fs_freq.append(0)
    for p, t in zip(ls_ls_presences, ls_ls_totals):
        if t > 0:
            ls_ls_freq.append(1 - p / t)
        else:
            ls_ls_freq.append(0)

    if interlaced:
        width = 0.8
        ls_all_labels = []
        palette = sns.color_palette().as_hex()
        palette_2 = sns.color_palette("husl", 8).as_hex()
        ls_colors_1 = []
        ls_colors_2 = []
        for i in range(len(ls_fs_freq)):
            if i == 0:
                ls_colors_1.append(palette[2])  # palette_2[2] or palette_2[3] or palette[8]
                # ls_colors_1.append('#efb100')
                ls_colors_2.append(palette[3])
                # ls_colors_2.append('#ef3d00')
            else:
                ls_colors_1.append(palette[9])
                ls_colors_2.append(palette[4])
            ls_all_labels.append(ls_fs_labels[i] + ' / ' + ls_ls_labels[i])

        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # ax.bar(x=ls_all_labels, height=ls_all_freq, width=0.8, align='center', color=ls_colors)
        ax.bar(x=np.arange(len(ls_fs_labels)) * 2 - width / 2, height=ls_fs_freq, width=width, align='center',
               color=ls_colors_1)
        ax.bar(x=np.arange(len(ls_ls_labels)) * 2 + width / 2, height=ls_ls_freq, width=width, align='center',
               color=ls_colors_2)
        ax.set_xticks(np.arange(len(ls_fs_labels)) * 2)
        ax.set_xticklabels(ls_all_labels)
        ax.set_xlabel('positions to the right/left of FS/LS')
        ax.set_ylabel('deletion frequencies')
        ax.set_title('CRISPRCasdb')
    else:
        palette = sns.color_palette().as_hex()
        palette = [palette[1] if (i == 0) else palette[0] for i in range(max_considered_distance + 1)]

        fig, ax = plt.subplots(1, 2, figsize=(20, 5))
        # fig.suptitle('FS/LS deletion frequencies')
        ax[0].bar(x=ls_fs_labels, height=ls_fs_freq, width=0.8, align='center', color=palette)
        ax[0].set_xlabel('positions to the right of FS')
        ax[0].set_ylabel('deletion frequencies')
        ax[0].set_title('FS+x deletion frequencies')

        ax[1].bar(x=ls_ls_labels[::-1], height=ls_ls_freq[::-1], width=0.8, align='center', color=palette[::-1])
        ax[1].set_xlabel('positions to the left of LS')
        ax[1].set_ylabel('deletion frequencies')
        ax[1].set_title('LS-x deletion frequencies')

        ymin = min(ax[0].get_ylim()[0], ax[1].get_ylim()[0])
        ymax = max(ax[0].get_ylim()[1], ax[1].get_ylim()[1])

        ax[0].set_ylim([ymin, ymax])
        ax[1].set_ylim([ymin, ymax])

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)

    sns.set_context("paper", font_scale=1)


def paper_plot_leaf_lengths(df, save_path, crispr_types=None):
    # sns.histplot(data=df, x='avg array length', binwidth=1, multiple='stack', hue='crispr type')
    df_plot = df.copy()
    array_count = sum([len(x) for x in df_plot['array names']])
    leaf_lengths = df['leaf array lengths']
    ls_leaf_lengths = []
    for i, d in enumerate(leaf_lengths):
        ls_leaf_lengths.append([int(l) for l in d.values()])
    df_plot['leaf array lengths'] = ls_leaf_lengths
    df_plot = df_plot.explode('leaf array lengths', ignore_index=True)
    nb_cols = 4
    if crispr_types is not None:
        df_plot = df_plot[df_plot['crispr type'].isin(crispr_types)]
        length_crispr_types = len(crispr_types)
    else:
        length_crispr_types = 1

    # print(df_plot)
    plt.figure(figsize=(15, 5))
    sns.histplot(data=df_plot, x='leaf array lengths', discrete=True)
    plt.xlabel(r'array length')
    plt.title(f'array length distribution [{array_count} arrays]')
    plt.tight_layout(pad=3.0)
    save_name = save_path + '_all.pdf' if crispr_types is not None else save_path
    plt.savefig(save_name)

    if crispr_types is not None:
        fig, ax = plt.subplots(length_crispr_types // nb_cols, nb_cols, figsize=(40, 15))
        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        for i, n in enumerate(crispr_types):
            df_data = df_plot[df_plot['crispr type'] == n]
            group_count = df_data.shape[0]
            row = i // nb_cols
            col = i % nb_cols
            sns.histplot(data=df_data, x='avg array length', discrete=True,
                         ax=ax[row, col])
            # sns.histplot(data=df, x='nb leafs (after combining)', discrete=True)
            plt.xlabel(r'Average array length per group')
            ax[row, col].set_title('CRISPR type ' + n + ' [%s groups]' % group_count)

        plt.tight_layout(pad=3.0)
        plt.savefig(save_path)
        plt.show()


def rename_columns_for_compatibility(df_results, inverse=False):
    rename_dict = {'ln_lh_idm': 'ln_lh_0', 'ln_lh_bdm': 'ln_lh_1', 'lh_idm': 'lh_0', 'lh_bdm': 'lh_1',
                   'deletion_rate_idm': 'loss_rate_0', 'deletion_rate_bdm': 'loss_rate_1',
                   'insertion_rate_idm': 'gain_rate_0', 'insertion_rate_bdm': 'gain_rate_1',
                   'test_statistic (-2*ln_lh_ratio)': '-2*ln_lh_ratio',
                   'alpha_bdm': 'alpha_1',
                   'per_spacer_deletion_rate_bdm': 'per spacer loss_rate_1',
                   'significance level': 'significance value',
                   'nb spacers in alignment': 'nb gained spacers',
                   'chi2_quantile': 'test_confidence',
                   'nb of leafs (before combining non-uniques)': 'nb leafs (before combining)',
                   'nb of leafs (after combining non-uniques)': 'nb leafs (after combining)',
                   'all max block deletion lengths': 'all max loss lengths',
                   'all max block deletion lengths (normalized)': 'all max loss lengths (normalized)',
                   'excluding young spacers deletion_rate_idm': 'ex_y_s loss_rate_0',
                   'excluding young spacers deletion_rate_bdm': 'ex_y_s loss_rate_1',
                   'excluding young spacers alpha_bdm': 'ex_y_s alpha_1',
                   'containing old spacers deletion_rate_idm': 'con_o_s loss_rate_0',
                   'containing old spacers deletion_rate_bdm': 'con_o_s loss_rate_1',
                   'containing old spacers alpha_bdm': 'con_o_s alpha_1',
                   'relative deletion positions': 'relative loss positions',
                   'reversed_ln_lh_idm': 'reversed_ln_lh_0', 'reversed_ln_lh_bdm': 'reversed_ln_lh_1',
                   'reversed_lh_idm': 'reversed_lh_0', 'reversed_lh_bdm': 'reversed_lh_1',
                   'lh_idm / reversed_lh_idm': 'lh_0 / reversed_lh_0',
                   'ln_lh_idm - reversed_ln_lh_idm': 'ln_lh_0 - reversed_ln_lh_0',
                   'lh_bdm / reversed_lh_bdm': 'lh_1 / reversed_lh_1',
                   'ln_lh_bdm - reversed_ln_lh_bdm': 'ln_lh_1 - reversed_ln_lh_1',
                   }
    if inverse:
        rename_dict = {v: k for k, v in rename_dict.items()}
    df_results = df_results.rename(columns=rename_dict)
    return df_results


if __name__ == '__main__':
    sns.set_theme()
    # data_folder = 'alex_data'
    # data_folder = os.path.join('0_result_folder', 'omer_data')
    # data_folder = os.path.join('0_result_folder', 'crispr_db')
    # data_folder = os.path.join('0_result_folder', 'sim_testing', 'test_data_from_pkl')
    data_folder = os.path.join('0_result_folder', '0_paper')
    # data_folder = 'real_data'
    # file_names = ['0_gr6_lr04_aNone', '1_gr6_lr04_a11', '2_gr6_lr04_a12', '3_gr6_lr04_a13', '4_gr6_lr04_a14', '5_gr6_lr04_a15']
    # file_names = ['0_rlp']
    # file_names = ['6_gr6_lr04_a10']
    # file_names = ['1_old_lh_full']
    # file_names = ['2_new_lh_full']
    # file_names = ['3_testing_tree_bl_plus_eps']
    # file_names = ['gb_crispr_type_new_lh_fct', 'gb_crispr_type_old_lh_fct']
    # file_names = ['testing_likelihood_trees_bl_old_lh', 'gb_crispr_type_old_lh_fct_fin_bl']
    # file_names = ['first_test_crispr_db']
    # file_names = ['determine_orientation_w_strand_pred']
    # file_names = ['pseudo_dna_fft']
    # file_names = ['determine_orientation_w_strand_pred_v2']
    # file_names = ['determine_orientation_w_strand_pred_v3']
    # file_names = ['determine_orientation_w_strand_pred_stdandfft']
    # file_names = ['determine_orientation_w_strand_pred_stdandfft_reindexed']
    # file_names = ['determine_orientation_w_strand_pred_stdandfft_reindexed_new_groups']
    # file_names = ['determine_orientation_w_stdandfft_reindexed_new_groups_combined_repeats_v2']
    # file_names = ['determine_orientation_w_stdandfft_reindexed_new_groups_combined_repeats_new_reindex']
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees']
    file_names = ['00_paper_real_data_optimal_estimation']
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees_alpha_bias']
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees_rho+alpha_bias_correction']
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees_rho+alpha_bias_correction_adjusted_test',
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees_rho+alpha_bias_correction_no_bias_test']
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees_alpha_bias_e-4']
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees_alpha_bias_e-5']
    # file_names = ['1_excluding_big_dataset_n_lh_fct_cg_trees_alpha_bias_e-7']
    for name in file_names:
        # path = os.path.join('pictures', 'real_data', name, '0_real_exp_protocol.csv')
        # path = os.path.join('pictures', data_folder, name, '0_real_exp_protocol.pkl')
        # path = os.path.join(data_folder, name, '0_protocol.pkl')
        # path = os.path.join(data_folder, name, '0_protocol_w_orientation.pkl')
        path = os.path.join(data_folder, name, '0_final_oriented_protocol.pkl')
        # path = os.path.join('pictures', 'alex_data', name, '0_real_exp_protocol.pkl')
        # path_boring = os.path.join('pictures', 'alex_data', name, '0_real_exp_protocol_boring.pkl')
        # save_path = os.path.join('pictures', 'real_data')
        # save_path = os.path.join('1_result_graphs', 'omer_data', name)
        # save_path = os.path.join('1_result_graphs', 'crispr_db_not_correct_orientation', name)
        # save_path = os.path.join('1_result_graphs', 'crispr_db', name)
        save_path = os.path.join('1_result_graphs', '0_paper', name)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)
        # df = pd.read_csv(path, index_col=0, converters={'relative loss positions': ast.literal_eval,
        #                                                 'all max loss lengths': ast.literal_eval})

        df = pd.read_pickle(path)

        # renaming in case of compatibility issues in column naming with old versions
        # df = rename_columns_for_compatibility(df)

        # df_boring = pd.read_pickle(path_boring)
        # df = pd.concat([df, df_boring])
        # print(df)
        sns.set_style("whitegrid")
        fig, ax = plt.subplots(figsize=(8, 6))
        sns.countplot(data=df, x='test result', order=[False, True])

        plt.xlabel('likelihood ratio test results')
        ax.set_xticklabels(['IDM', 'BDM'])

        plt.tight_layout(pad=1.0)
        plt.savefig(os.path.join(save_path, '_'.join([name, 'test_results_count.pdf'])))
        plt.show()

        plot_paper_test_results_vs_type(df,
                                        os.path.join(save_path,
                                                     '_'.join([name,
                                                               'paper_test_results_vs_cas_type.pdf'])),
                                        order=None, x='crispr type', rotate_90=False,
                                        count_threshold=10,
                                        figsize=(15, 5))

        paper_plot_leaf_lengths(df, os.path.join(save_path, '_'.join([name, 'leaf_lengths.pdf'])), )

        plot_for_paper_estimated_per_spacer_loss_rate(df, os.path.join(save_path,
                                                                       '_'.join([name,
                                                                                 'paper_per_spacer_loss_rate.pdf'])),
                                                      order=None, x='crispr type', rotate_90=False, count_threshold=10)

        plot_paper_fs_ls_loss_freq(df, os.path.join(save_path, '_'.join([name,
                                                                         'paper_fs_ls_loss_freq_interlaced.pdf'])),
                                   max_considered_distance=5, binning=False, interlaced=True)

        plot_paper_fs_ls_loss_freq(df, os.path.join(save_path, '_'.join([name, 'paper_fs_ls_loss_freq.pdf'])),
                                   max_considered_distance=5, binning=False)
        plot_paper_fs_ls_loss_freq(df, os.path.join(save_path, '_'.join([name, 'paper_fs_ls_loss_freq_binned.pdf'])),
                                   max_considered_distance=5, binning=True)

        plot_poster_array_orientation_agreement(df,
                                                os.path.join(save_path, '_'.join([name,
                                                                                  'poster_orientation_solo_agreement_arrays.pdf'])),
                                                filter_crispr_types=['CAS-TypeI-B', 'CAS-TypeI-C', 'CAS-TypeI-E',
                                                                     'CAS-TypeI-F',
                                                                     'CAS-TypeII-A', 'CAS-TypeII-C'])

        plot_for_poster_estimated_per_spacer_loss_rate_lrt(df,
                                                           os.path.join(save_path,
                                                                        '_'.join([name,
                                                                                  'poster_per_spacer_loss_rate.pdf'])),
                                                           order=None, x='crispr type', rotate_90=False,
                                                           count_threshold=10)
        plot_for_poster_estimated_parameters_vs_genus(df,
                                                      os.path.join(save_path,
                                                                   '_'.join([name,
                                                                             'poster_estimated_parameters_vs_genus.pdf'])),
                                                      order=None,
                                                      x='group species',
                                                      rotate_90=True,
                                                      count_threshold=10)
        plot_estimated_parameters(df,
                                  os.path.join(save_path,
                                               '_'.join([name, 'estimated_para_group_species_count_threshold.pdf'])),
                                  order=None,
                                  x='group species',
                                  rotate_90=True,
                                  figsize=(50, 5),
                                  count_threshold=10,
                                  log_transform=False,
                                  remove_extreme_outliers=True)

        plot_rel_pos(df, os.path.join(save_path, '_'.join([name, 'rel_pos_hist.pdf'])))

        plot_can_predict(df, os.path.join(save_path, '_'.join([name, 'poster_orientation_arrays.pdf'])),
                         filter_crispr_types=['CAS-TypeI-B', 'CAS-TypeI-C', 'CAS-TypeI-E', 'CAS-TypeI-F',
                                              'CAS-TypeII-A', 'CAS-TypeII-C']
                         )
        plot_can_predict(df, os.path.join(save_path, '_'.join([name, 'poster_orientation_groups.pdf'])), by_group=True,
                         filter_crispr_types=['CAS-TypeI-B', 'CAS-TypeI-C', 'CAS-TypeI-E', 'CAS-TypeI-F',
                                              'CAS-TypeII-A', 'CAS-TypeII-C'])

        # plot_orientation_likelihood_distribution(df, os.path.join(save_path, '_'.join([name,
        #                                                                                'orientation_lh_dist.pdf'])),
        #                                          log=False)
        plot_rel_pos_freq(df, os.path.join(save_path, '_'.join([name, 'rel_loss_frequencies.pdf'])))

        plot_fes_les_loss_freq(df,
                               os.path.join(save_path, '_'.join([name, 'fes_les_presence_freq.pdf'])))
        plot_fes_les_loss_freq(df,
                               os.path.join(save_path, '_'.join([name, 'global_fes_les_presence_freq.pdf'])),
                               use_global=True)

        plot_orientation_likelihood_distribution(df, os.path.join(save_path, '_'.join([name,
                                                                                       'orientation_ln_lh_dist.pdf'])),
                                                 log=True)
        plot_orientation_likelihood_distribution(df, os.path.join(save_path,
                                                                  '_'.join([name,
                                                                            'orientation_ln_lh_dist_wo_3.pdf'])),
                                                 log=True,
                                                 filter_3=True)

        plot_max_loss_lengths(df, os.path.join(save_path, '_'.join([name, 'max_loss_length.pdf'])))
        plot_max_loss_lengths_norm(df, os.path.join(save_path, '_'.join([name, 'max_loss_length_norm.pdf'])))

        plot_count_nb_unique_leafs(df, os.path.join(save_path, '_'.join([name, 'nb_unique_leafs.pdf'])))
        plot_count_nb_leafs(df, os.path.join(save_path, '_'.join([name, 'nb_leafs.pdf'])))
        plot_avg_array_length(df, os.path.join(save_path, '_'.join([name, 'avg_array_length.pdf'])))
        plot_nb_unique_spacers(df, os.path.join(save_path, '_'.join([name, 'nb_gained_spacers.pdf'])))

        plot_type_counts(df, os.path.join(save_path, '_'.join([name, 'crispr_type_counts.pdf'])))
        plot_type_counts(df, os.path.join(save_path, '_'.join([name, 'crispr_type_counts_nb_arrays.pdf'])),
                         nb_arrays=True)
        plot_estimated_parameters(df, os.path.join(save_path, '_'.join([name, 'estimated_para.pdf'])),
                                  # order=['CAS-I-C', 'CAS-I-E', 'CAS-I-F', 'CAS-II-A', 'CAS-II-C', 'CAS-III-A'],
                                  order=None,
                                  x='crispr type',
                                  rotate_90=True
                                  )
        plot_estimated_parameters(df, os.path.join(save_path, '_'.join([name, 'estimated_para_ch_pl.pdf'])),
                                  order=None,
                                  x='chromosome_plasmid',
                                  )
        plot_group_species_counts(df, os.path.join(save_path, '_'.join([name, 'group_species_counts.pdf'])))
        plot_group_kingdom_counts(df, os.path.join(save_path, '_'.join([name, 'group_kingdom_counts.pdf'])))
        # plot_individual_species_counts(df, os.path.join(save_path, '_'.join([name, 'ind_species_counts.pdf'])))
        plot_individual_genus_counts(df, os.path.join(save_path, '_'.join([name, 'ind_genus_counts.pdf'])))
        plot_estimated_parameters(df, os.path.join(save_path, '_'.join([name, 'estimated_para_group_species.pdf'])),
                                  order=None,
                                  x='group species',
                                  rotate_90=True,
                                  figsize=(50, 5),
                                  log_transform=False)
        plot_for_paper_estimated_gain_rate(df, os.path.join(save_path, '_'.join([name,
                                                                                 'paper_estimated_gain_rate.pdf'])),
                                           order=None, x='crispr type', rotate_90=False, count_threshold=10)
        plot_estimated_parameters(df, os.path.join(save_path, '_'.join([name, 'estimated_para_count_threshold.pdf'])),
                                  # order=['CAS-I-C', 'CAS-I-E', 'CAS-I-F', 'CAS-II-A', 'CAS-II-C', 'CAS-III-A'],
                                  order=None,
                                  x='crispr type',
                                  count_threshold=10,
                                  rotate_90=False,
                                  )
        plot_estimated_parameters(df, os.path.join(save_path, '_'.join(
            [name, 'estimated_para_count_threshold_atleast4_strains.pdf'])),
                                  # order=['CAS-I-C', 'CAS-I-E', 'CAS-I-F', 'CAS-II-A', 'CAS-II-C', 'CAS-III-A'],
                                  order=None,
                                  x='crispr type',
                                  count_threshold=10,
                                  rotate_90=False,
                                  exclude_less_strains=4,
                                  )
        plot_count_orientation(df,
                               os.path.join(save_path, '_'.join([name, 'orientation_count.pdf'])))

        plot_orientation_stats(df,
                               os.path.join(save_path, '_'.join([name, 'orientation_stats.pdf'])))
        plot_orientation_stats(df,
                               os.path.join(save_path, '_'.join([name, 'orientation_stats_groups.pdf'])),
                               by_group=True)

        large_types = ['CAS-TypeI-B', 'CAS-TypeI-C', 'CAS-TypeI-E', 'CAS-TypeI-F',
                       'CAS-TypeII-A', 'CAS-TypeII-C',
                       'CAS-TypeIII-A']
        # plot_max_loss_lengths(df, os.path.join(save_path, '_'.join([name, 'max_loss_length_large_types.pdf'])),
        #                       crispr_types=large_types)

        # plot_max_loss_lengths_norm(df,
        #                            os.path.join(save_path, '_'.join([name, 'max_loss_length_norm_large_types.pdf'])),
        #                            crispr_types=large_types)
        plot_count_nb_unique_leafs(df,
                                   os.path.join(save_path, '_'.join([name, 'nb_unique_leafs_large_types.pdf'])),
                                   crispr_types=large_types)
        plot_count_nb_leafs(df, os.path.join(save_path, '_'.join([name, 'nb_leafs_large_types.pdf'])),
                            crispr_types=large_types)
        plot_avg_array_length(df, os.path.join(save_path, '_'.join([name, 'avg_array_length_large_types.pdf'])),
                              crispr_types=large_types)

        plot_orientation_stats(df,
                               os.path.join(save_path, '_'.join([name, 'orientation_stats_large_types.pdf'])),
                               filter_crispr_types=large_types)
        plot_orientation_likelihood_distribution(df,
                                                 os.path.join(save_path,
                                                              '_'.join([name,
                                                                        'orientation_ln_lh_dist_large_types.pdf'])),
                                                 filter_crispr_types=large_types, log=True)
        plot_orientation_stats(df,
                               os.path.join(save_path, '_'.join([name, 'orientation_stats_groups_large_types.pdf'])),
                               by_group=True, filter_crispr_types=large_types)
        plot_orientation_likelihood_distribution(df,
                                                 os.path.join(save_path,
                                                              '_'.join([name,
                                                                        'orientation_ln_lh_dist_large_types.pdf'])),
                                                 filter_crispr_types=large_types, log=True)
        plot_orientation_stats(df,
                               os.path.join(save_path, '_'.join([name, 'orientation_stats_wo_3.pdf'])),
                               filter_3=True)
        plot_orientation_stats(df,
                               os.path.join(save_path, '_'.join([name, 'orientation_stats_groups_mt3arr.pdf'])),
                               filter_more_than_3_unique=True, by_group=True)

        for t in large_types:
            plot_orientation_likelihood_distribution(df,
                                                     os.path.join(save_path,
                                                                  '_'.join([name,
                                                                            'orientation_ln_lh_dist_' + t + '.pdf'])),
                                                     filter_crispr_types=[t], log=True)
            plot_orientation_stats(df,
                                   os.path.join(save_path, '_'.join([name, 'orientation_stats_groups_' + t + '.pdf'])),
                                   by_group=True, filter_crispr_types=[t])
            # plot_fes_les_loss_freq(df,
            #                        os.path.join(save_path, '_'.join([name, 'fes_les_presence_freq_' + t + '.pdf'])),
            #                        crispr_type=t)
            # plot_fes_les_loss_freq(df,
            #                        os.path.join(save_path, '_'.join([name,
            #                                                          'global_fes_les_presence_freq_' + t + '.pdf'])),
            #                        use_global=True, crispr_type=t)
            plot_rel_pos_freq(df, os.path.join(save_path, '_'.join([name, 'rel_loss_frequencies_' + t + '.pdf'])),
                              crispr_type=t)
            plot_rel_pos(df, os.path.join(save_path, '_'.join([name, 'rel_pos_hist_' + t + '.pdf'])),
                         crispr_type=t)
