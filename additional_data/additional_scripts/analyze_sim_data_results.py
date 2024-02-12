import collections
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import ast
import numpy as np
import copy

BIN_WIDTH = 0.02


def read_csv_string_ls_to_ls(path, column):
    return pd.read_csv(path, converters={column: ast.literal_eval})


def plot_orientation_threshold_roc_curve(df, ls_threshold_values, save_path):
    # filter?
    df_plot = copy.deepcopy(df)

    def decision_fct(ln_lh_diff, decision_boundary=0):
        decision = 'nd'
        if ln_lh_diff < -decision_boundary:
            decision = True
        elif ln_lh_diff > decision_boundary:
            decision = False
        return decision

    df_fct_data = pd.DataFrame(columns=['threshold', 'decisiveness', 'precision'])
    df_fct_data['threshold'] = ls_threshold_values

    ls_decisiveness, ls_precision = [], []

    for threshold in ls_threshold_values:
        df_plot[str(threshold)] = [decision_fct(v, decision_boundary=threshold)
                                   for v in df_plot['ln_lh_1 - reversed_ln_lh_1']]
        ls_decisiveness.append(1 - df_plot[str(threshold)].value_counts(normalize=True)['nd'])
        # precision is only the "False" predictions, since all simulated data should be oriented forward,
        # i.e. correctly, i.e. no reversion is recommended.
        print(threshold, df_plot[str(threshold)].value_counts())
        df_prec = df_plot[df_plot[str(threshold)] != 'nd']
        prec_value_count = df_prec[str(threshold)].value_counts(normalize=True)[False] if not df_prec.empty else 0
        print(df_prec[str(threshold)].value_counts(normalize=True))
        ls_precision.append(prec_value_count)

    df_fct_data['decisiveness'] = ls_decisiveness
    df_fct_data['precision'] = ls_precision
    print(df_fct_data)

    plt.figure()
    sns.lineplot(x='threshold', y='decisiveness', data=df_fct_data, label='decisiveness')
    sns.lineplot(x='threshold', y='precision', data=df_fct_data, label='precision')
    plt.xlabel('Orientation threshold')
    plt.ylabel('Proportion')
    plt.ylim(0, 1)
    plt.legend()
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


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


def plot_rel_pos_freq(df, save_path, nb_cols=3, bin_width=0.01, sim=None, ):
    # array_length_bins = [0, 10, 15, 20, 30, np.infty]
    # names = ['<10', '10-15', '15-20', '20-30', '30+']
    # bin_width = [7, 13, 17, 25, 31, 31]
    array_length_bins = [0, 10, 15, 25, 30, np.infty]
    names = ['<10', '10-15', '15-25', '25-30', '30+']
    bin_width = [7, 13, 15, 25, 31, 31]
    df['avg array length bins'] = pd.cut(df['avg array length'], array_length_bins, labels=names)
    df_start = df.copy()
    names.append('all')
    dict_df_plot = {}
    dict_nb_groups_per_bin = {'all': df_start.shape[0]}
    # print(df_start)
    for j, array_length_label in enumerate(names):
        if array_length_label == 'all':
            df_plot = df_start
        else:
            df_plot = df_start[df_start['avg array length bins'] == array_length_label]
        positions = []
        if sim == 'branch_based_sim':
            ex_spacers_key = 'sim branch based nb of existent spacers'
        else:
            ex_spacers_key = 'sim nb of existent spacers' if sim else 'nb of existent spacers'
        for ls_ex_spacers in df_plot[ex_spacers_key].values:
            # print(ls_ex_spacers)
            for nb in ls_ex_spacers:
                if nb > 1:
                    positions += [i / (nb - 1) for i in range(nb)]
                else:
                    positions += [1.0]
        dict_pos_counts = collections.Counter(positions)
        all_rel_loss_pos = []
        if sim == 'branch_based_sim':
            rlp_key = 'sim. branch based relative loss positions'
            new_ls = []
            for idx in df_plot.index:
                new_ls.append([v for val in df_plot.loc[idx, rlp_key] for v in val])
            df_plot[rlp_key] = new_ls
        else:
            rlp_key = 'sim. relative loss positions' if sim else 'relative loss positions'
        for ls_rel_loss_pos in df_plot[rlp_key].values:
            all_rel_loss_pos += ls_rel_loss_pos
        # handling empty case...
        if not all_rel_loss_pos:
            df_freq = pd.DataFrame()
            dict_df_plot[array_length_label] = df_freq
            dict_nb_groups_per_bin[array_length_label] = df_plot.shape[0]
            continue
        dict_rel_loss_counts = collections.Counter(all_rel_loss_pos)
        bins = list(np.linspace(0, 1, num=bin_width[j]))
        p_pos_counts_bins, bin_labels = digitize(list(dict_pos_counts.keys()), bins)
        p_rel_loss_counts_bins, _ = digitize(list(dict_rel_loss_counts.keys()), bins)
        bins = bin_labels
        binned_dict_pos_counts = {}
        binned_rel_loss_counts = {}
        for b_idx, count in zip(p_pos_counts_bins, dict_pos_counts.values()):
            binned_dict_pos_counts[bins[b_idx]] = binned_dict_pos_counts.get(bins[b_idx], 0) + count
        for b_idx, count in zip(p_rel_loss_counts_bins, dict_rel_loss_counts.values()):
            binned_rel_loss_counts[bins[b_idx]] = binned_rel_loss_counts.get(bins[b_idx], 0) + count

        dict_freq = {p: binned_rel_loss_counts.get(p, 0) / p_count for p, p_count in binned_dict_pos_counts.items()}
        dict_freq = [{'positions': p, 'relative loss frequencies': freq} for p, freq in dict_freq.items()]
        # print('dict_freq', dict_freq)
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
            plt.xlabel('relative array position' if sim else 'relative array position')
            plt.ylabel('deletion frequency' if sim else 'deletion frequency')
            plt.xticks(rotation=90)
            plt.title('Avg. array length bin ' + n + ' [%s groups]' % dict_nb_groups_per_bin[n])

            # plt.title('Deletion frequencies of groups with Avg. array length bin ' + n + '[%s groups]'
            #           % dict_nb_groups_per_bin[n])
            plt.tight_layout(pad=1.0)
            plt.savefig(save_path + '15-25_sim.pdf' if sim else save_path + '15-25.pdf')

        sns.barplot(data=df_freq, x='positions', y='relative loss frequencies', ax=ax[row, col], palette=palette)
        # ax[row, col].bar(x=df_freq['positions'], height=df_freq['relative loss frequencies'], width=bin_width,
        #                  align='center')

        ax[row, col].tick_params(axis='x', rotation=90)
        ax[row, col].set_xlabel('relative array position' if sim else 'relative array position')
        ax[row, col].set_ylabel('deletion frequency' if sim else 'deletion frequency')
        # sns.histplot(data=df_plot[df_plot['avg array length bins'] == n], x='relative loss positions',
        #              binwidth=BIN_WIDTH, ax=ax[row, col])
        ax[row, col].set_title('Avg. array length bin ' + n + ' [%s groups]' % dict_nb_groups_per_bin[n])
    nb_of_plots = int(np.ceil(len(names) / nb_cols)) * nb_cols

    for j in range(nb_of_plots - len(names)):
        row = - j // nb_cols - 1
        col = - j % nb_cols - 1
        ax[row, col].axis('off')

    fig.tight_layout(pad=3.0)
    fig.savefig(save_path)
    fig.show()


def plot_rel_pos(df, save_path, nb_cols=3):
    df_plot = df.explode('relative loss positions', ignore_index=True)
    bins = [0, 10, 15, 20, 30, np.infty]
    names = ['<10', '10-15', '15-20', '20-30', '30+']
    df_plot['avg array length bins'] = pd.cut(df_plot['avg array length'], bins, labels=names)
    # print('0', df_plot[(df_plot['relative loss positions'] == 0)].shape[0])
    # print('0.001', df_plot[(df_plot['relative loss positions'] > 0) & (df_plot['relative loss positions'] <= 0.001)].shape[0])
    # print(df_plot)
    # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    fig, ax = plt.subplots(len(names) // nb_cols + 1, nb_cols, figsize=(30, 15))
    for i, n in enumerate(names):
        row = i // nb_cols
        col = i % nb_cols

        sns.histplot(data=df_plot[df_plot['avg array length bins'] == n], x='relative loss positions',
                     binwidth=BIN_WIDTH, ax=ax[row, col])
        ax[row, col].set_title('Avg. array length bin ' + n)
    # sns.histplot(data=df_plot, x='relative loss positions', binwidth=BIN_WIDTH, multiple='stack',
    #              hue='avg array length bins')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_sim_rel_pos(df, save_path, nb_cols=3):
    df_plot = df.explode('sim. relative loss positions', ignore_index=True)
    bins = [0, 10, 15, 20, 30, np.infty]
    bin_width = [7, 13, 17, 25, 31, 31]
    names = ['<10', '10-15', '15-20', '20-30', '30+']
    df_plot['avg array length bins'] = pd.cut(df_plot['avg array length'], bins, labels=names)
    # print('0', df_plot[(df_plot['sim. relative loss positions'] == 0)].shape[0])
    # print('0.001', df_plot[(df_plot['sim. relative loss positions'] > 0) & (df_plot['sim. relative loss positions'] <= 0.001)].shape[0])
    # print(df_plot)
    # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    names.append('all')
    dict_nb_groups_per_bin = {'all': df_plot.shape[0]}

    fig, ax = plt.subplots(len(names) // nb_cols + 1, nb_cols, figsize=(30, 15))
    for i, n in enumerate(names):
        if n == 'all':
            df_plot = df_plot
        else:
            df_plot = df_plot[df_plot['avg array length bins'] == n]
        # palette = ['tab:orange' if (i == 0 or i == len(df_freq['positions']) - 1) else 'tab:blue'
        #            for i, _ in enumerate(df_freq['positions'])]

        row = i // nb_cols
        col = i % nb_cols

        sns.histplot(data=df_plot, x='sim. relative loss positions',
                     binwidth=bin_width[i], ax=ax[row, col])
        ax[row, col].set_title('Avg. array length bin ' + n)
    # sns.histplot(data=df_plot, x='sim. relative loss positions', binwidth=BIN_WIDTH, multiple='stack',
    #              hue='avg array length bins')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_max_loss_lengths(df, save_path, select_frag_par=None):
    df_plot = df.copy()
    fig = plt.figure(figsize=(10, 5))
    if select_frag_par is not None:
        # df = df.loc[df['frag. loss param.'] == select_frag_par]
        df_plot = df_plot.loc[df['loss model param.'] == select_frag_par]
    df_plot = df_plot.explode('all max loss lengths', ignore_index=True)
    # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    sns.histplot(data=df_plot, x='all max loss lengths', stat='proportion', discrete=True)
    plt.xlabel('adjacent deletion lengths')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_max_loss_lengths_norm(df, save_path, select_frag_par=None):
    df_plot = df.copy()
    fig = plt.figure(figsize=(10, 5))
    if select_frag_par is not None:
        # df = df.loc[df['frag. loss param.'] == select_frag_par]
        df_plot = df_plot.loc[df['loss model param.'] == select_frag_par]
    df_plot = df_plot.explode('all max loss lengths (normalized)', ignore_index=True)
    # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    sns.histplot(data=df_plot, x='all max loss lengths (normalized)', stat='probability')
    plt.xlabel('adjacent deletion lengths (normalized)')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_test_results_vs_frag_para(df, save_path):
    df_1 = df.groupby('loss model param.')['test result'].value_counts(normalize=True)
    df_1 = df_1.rename('test result normalized').reset_index()
    sns.catplot(data=df_1, x='loss model param.', y='test result normalized', kind='bar', hue='test result')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_frag_para_bin_plot(df, save_path):
    sns.displot(data=df, x='loss model param.', binwidth=2, hue='test result', multiple='stack')
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_nb_gained_spacers(df, save_path):
    sns.histplot(data=df, x='nb gained spacers', discrete=True)
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_alternatives_estimated_parameters(df, save_path, showfliers=False):
    df_plot = df[df['all max loss lengths'].astype(bool)]
    fig, ax = plt.subplots(3, 3, figsize=(20, 15))
    unique_sim_lr = pd.unique(df_plot['sim. loss rate'])
    # sns.countplot(data=df, x='loss model param.', hue='test result', ax=ax[0], )
    # # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    # ax[0].set(title='Test results')
    # ax[0].legend(title='Test result')
    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_0', hue=None, ax=ax[0, 0], showfliers=showfliers,
                )
    ax[0, 0].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red')
    ax[0, 0].set(ylabel=r'loss rate $\hat{\rho}_0$', title=r'Independent loss model ($\alpha = 1$)')

    # ax[1].legend(title='Test result')
    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_1', hue=None, ax=ax[1, 0], showfliers=showfliers,
                )
    # sns.stripplot(data=df, x='loss model param.', y='sim. loss rate', ax=ax[2], hue=None, marker='_', facecolor='red',
    #               edgecolor=None, size=100)
    for i, val in enumerate(unique_sim_lr):
        ax[1, 0].axhline(val, i * 1 / len(unique_sim_lr), (i + 1) * 1 / len(unique_sim_lr), color='red')
    ax[1, 0].set(ylabel=r'loss rate $\hat{\rho}_1$')  # , title=r'Block loss model: initialization rate'
    # ax[2].legend(title='Test result')
    sns.boxplot(data=df_plot, x='loss model param.', y='alpha_1', hue=None, ax=ax[2, 0], showfliers=showfliers, )
    unique_sim_alpha = pd.unique(df_plot['loss model param.'].dropna())
    for i, val in enumerate(unique_sim_alpha):
        ax[2, 0].axhline(val, i * 1 / len(unique_sim_alpha), (i + 1) * 1 / len(unique_sim_alpha), color='red')
    # unique_sim_alpha = pd.unique(df['loss model param.'])
    # for i, val in enumerate(unique_sim_alpha):
    #     ax[3].axhline(val, i*1/len(unique_sim_alpha), (i + 1)*1/len(unique_sim_alpha), color='red')
    ax[2, 0].set(ylabel=r'expected block length $\hat{\alpha}$')  # , title=r'Block loss model: expected block length'
    # ax[3].legend(title='Test result')
    # plt.figtext(0.76, 0.95, 'Block loss model', ha='center', va='center')

    sns.boxplot(data=df_plot, x='loss model param.', y='ex_y_s loss_rate_0', hue=None, ax=ax[0, 1],
                showfliers=showfliers,
                )
    ax[0, 1].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red')

    ax[0, 1].set(ylabel=r'ex ys del $\hat{\rho}_0$', title=r'excluding young spacer deletions ($\alpha = 1$)')
    # ax[1].legend(title='Test result')
    sns.boxplot(data=df, x='loss model param.', y='ex_y_s loss_rate_1', hue=None, ax=ax[1, 1], showfliers=showfliers,
                )
    # sns.stripplot(data=df, x='loss model param.', y='sim. loss rate', ax=ax[2], hue=None, marker='_', facecolor='red',
    #               edgecolor=None, size=100)
    unique_sim_lr = pd.unique(df['sim. loss rate'])
    for i, val in enumerate(unique_sim_lr):
        ax[1, 1].axhline(val, i * 1 / len(unique_sim_lr), (i + 1) * 1 / len(unique_sim_lr), color='red')
    ax[1, 1].set(ylabel=r'ex ys del $\hat{\rho}_1$')  # , title=r'Block loss model: initialization rate'
    # ax[2].legend(title='Test result')
    sns.boxplot(data=df, x='loss model param.', y='ex_y_s alpha_1', hue=None, ax=ax[2, 1], showfliers=showfliers, )

    for i, val in enumerate(unique_sim_alpha):
        ax[2, 1].axhline(val, i * 1 / len(unique_sim_alpha), (i + 1) * 1 / len(unique_sim_alpha), color='red')
    # unique_sim_alpha = pd.unique(df['loss model param.'])
    # for i, val in enumerate(unique_sim_alpha):
    #     ax[3].axhline(val, i*1/len(unique_sim_alpha), (i + 1)*1/len(unique_sim_alpha), color='red')
    ax[2, 1].set(ylabel=r'ex ys del $\hat{\alpha}$')  # , title=r'Block loss model: expected block length'
    sns.boxplot(data=df, x='loss model param.', y='con_o_s loss_rate_0', hue=None, ax=ax[0, 2], showfliers=showfliers,
                )
    ax[0, 2].set(ylabel=r'con os $\hat{\rho}_0$', title=r'Contains old spacer ($\alpha = 1$)')
    ax[0, 2].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red')
    # ax[1].legend(title='Test result')
    sns.boxplot(data=df, x='loss model param.', y='con_o_s loss_rate_1', hue=None, ax=ax[1, 2], showfliers=showfliers,
                )
    # sns.stripplot(data=df, x='loss model param.', y='sim. loss rate', ax=ax[2], hue=None, marker='_', facecolor='red',
    #               edgecolor=None, size=100)
    unique_sim_lr = pd.unique(df['sim. loss rate'])
    for i, val in enumerate(unique_sim_lr):
        ax[1, 2].axhline(val, i * 1 / len(unique_sim_lr), (i + 1) * 1 / len(unique_sim_lr), color='red')
    ax[1, 2].set(ylabel=r'con os $\hat{\rho}_1$')  # , title=r'Block loss model: initialization rate'
    # ax[2].legend(title='Test result')
    sns.boxplot(data=df, x='loss model param.', y='con_o_s alpha_1', hue=None, ax=ax[2, 2], showfliers=showfliers, )

    for i, val in enumerate(unique_sim_alpha):
        ax[2, 2].axhline(val, i * 1 / len(unique_sim_alpha), (i + 1) * 1 / len(unique_sim_alpha), color='red')
    # unique_sim_alpha = pd.unique(df['loss model param.'])
    # for i, val in enumerate(unique_sim_alpha):
    #     ax[3].axhline(val, i*1/len(unique_sim_alpha), (i + 1)*1/len(unique_sim_alpha), color='red')
    ax[2, 2].set(ylabel=r'con os $\hat{\alpha}$')  # , title=r'Block loss model: expected block length'
    fig.tight_layout(pad=1.0)
    for i in range(ax.shape[0]):
        ylim_upper = 0
        ylim_lower = 0
        for j in range(ax.shape[1]):
            ylim_upper = max(ylim_upper, ax[i, j].get_ylim()[1])
            ylim_lower = max(ylim_lower, ax[i, j].get_ylim()[0])
        for j in range(ax.shape[1]):
            ax[i, j].set_ylim([ylim_lower, ylim_upper])
    plt.show()
    fig.savefig(save_path)


def paper_plot_estimated_parameters(df, save_path, showfliers=True):
    df_plot = df.copy()
    ls_alpha = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,
                3, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
                4, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9,
                5]
    sns.set_context("paper", font_scale=3.5)
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]
    fig = plt.figure(figsize=(40, 15))
    gs = GridSpec(2, 10, figure=fig)
    # gs.update(hspace=3)
    # title = fig.add_sub9dplot(gs[0, :])
    ax = [fig.add_subplot(gs[0, :2]), fig.add_subplot(gs[0, 2:6]), fig.add_subplot(gs[0, 6:]),
          fig.add_subplot(gs[1, :5]), fig.add_subplot(gs[1, 5:])]

    df_plot = df_plot[df_plot['loss model param.'].isin(ls_alpha)]
    unique_sim_lr = pd.unique(df_plot['sim. loss rate'])

    df_count_plot = df_plot[df_plot['loss model param.'].isin(ls_alpha[:21])]
    # sns.countplot(data=df_count_plot, x='loss model param.', hue='test result', ax=ax[0], palette='colorblind')
    df_count_plot['loss model param.'] = df_count_plot['loss model param.'].astype(str)
    x, y = 'loss model param.', 'test result'

    df_count_plot[y] = df_count_plot[y].replace({True: 'BDM', False: 'IDM'})
    df_count_plot = df_count_plot.groupby(x)[y].value_counts(normalize=True)
    df_count_plot = df_count_plot.mul(100)
    df_count_plot = df_count_plot.rename('percent').reset_index()
    sns.lineplot(data=df_count_plot, x=x, y='percent', hue=y, ax=ax[0], linewidth=4)
    # sns.histplot(data=df_count_plot, x='loss model param.', hue='test result', ax=ax[0], multiple='dodge',
    #              stat='probability', common_norm=False)
    # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]

    # ax[0].set(title='Test results')
    h, l = ax[0].get_legend_handles_labels()
    ax[0].legend(h, l, title='test result')
    ax[0].set(xlabel=r'expected block length $\alpha$', )  # of sim.
    for legobj in ax[0].get_legend().legend_handles:
        legobj.set_linewidth(4.0)

    for ind, label in enumerate(ax[0].get_xticklabels()):
        if ind % 4 == 0:
            label.set_visible(True)
        else:
            label.set_visible(False)

    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_0', hue=None, ax=ax[1], showfliers=showfliers, )
    ax[1].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red', linewidth=4)

    ax[1].set(xlabel=r'expected block length $\alpha$',  # of simulation
              ylabel='',
              title=r'estimated IDM deletion rate $\hat{\rho}_I$',
              # title=r'Independent Deletion Model ($\alpha = 1$)'
              )

    for ind, label in enumerate(ax[1].get_xticklabels()):
        if ind % 4 == 0:
            label.set_visible(True)
        else:
            label.set_visible(False)
    # ax[1].legend(title='Test result')
    # sns.stripplot(data=df, x='loss model param.', y='sim. loss rate', ax=ax[2], hue=None, marker='_', facecolor='red',
    #               edgecolor=None, size=100)

    df_plot.loc[:, 'loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_1_times_alpha_1',
                hue=None, ax=ax[2], showfliers=showfliers, )
    ax[2].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red', linewidth=4)
    ax[2].set(xlabel=r'expected block length $\alpha$',  # of simulation
              ylabel='',
              title=r'estimated BDM per spacer deletion rate $\hat{\rho}_B \cdot \hat{\alpha}$')

    for ind, label in enumerate(ax[2].get_xticklabels()):
        if ind % 4 == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)

    ylim_upper = max(ax[1].get_ylim()[1], ax[2].get_ylim()[1])
    ylim_lower = min(ax[1].get_ylim()[0], ax[2].get_ylim()[0])
    ax[1].set_ylim([ylim_lower, ylim_upper])
    ax[2].set_ylim([ylim_lower, ylim_upper])

    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_1', hue=None, ax=ax[3], showfliers=showfliers,
                )
    for i, val in enumerate(unique_sim_lr):
        ax[3].axhline(val, i * 1 / len(unique_sim_lr), (i + 1) * 1 / len(unique_sim_lr), color='red', linewidth=4)
    ax[3].set(xlabel=r'expected block length $\alpha$',  # of simulation
              ylabel='',
              title=r'estimated BDM deletion rate $\hat{\rho}_B$')  # , title=r'Block loss model: initialization rate'
    # ax[2].legend(labels=[r'simulation $\rho_B$'], )
    # ax[2].legend(title='Test result')

    for ind, label in enumerate(ax[3].get_xticklabels()):
        if ind % 2 == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)

    sns.boxplot(data=df_plot, x='loss model param.', y='alpha_1', hue=None, ax=ax[4], showfliers=showfliers, )

    unique_sim_alpha = pd.unique(df_plot['loss model param.'].dropna())
    for i, val in enumerate(unique_sim_alpha):
        ax[4].axhline(val, i * 1 / len(unique_sim_alpha), (i + 1) * 1 / len(unique_sim_alpha), color='red', linewidth=4)
    # unique_sim_alpha = pd.unique(df['loss model param.'])
    # for i, val in enumerate(unique_sim_alpha):
    #     ax[3].axhline(val, i*1/len(unique_sim_alpha), (i + 1)*1/len(unique_sim_alpha), color='red')
    ax[4].set(xlabel=r'expected block length $\alpha$',  # of simulation
              ylabel='',
              title=r'estimated BDM expected block length $\hat{\alpha}$')
    # ax[3].legend(title='Test result')

    for ind, label in enumerate(ax[4].get_xticklabels()):
        if ind % 2 == 0:  # every 10th label is kept
            label.set_visible(True)
        else:
            label.set_visible(False)
    # plt.figtext(0.76, 0.95, 'Block Deletion Model', ha='center', va='center')
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)
    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1)


def plot_estimated_parameters(df, save_path, showfliers=True):
    """
    :param df:
    :param save_path:
    :param showfliers:
    :return:
    """
    df_plot = df[df['all max loss lengths'].astype(bool)]
    unique_sim_lr = pd.unique(df_plot['sim. loss rate'])
    fig, ax = plt.subplots(1, 6, figsize=(35, 5))
    sns.countplot(data=df_plot, x='loss model param.', hue='test result', ax=ax[0], )
    # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]
    ax[0].set(title='Test results')
    ax[0].legend(title='Test result', labels=['IDM', 'BDM'])
    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_0', hue=None, ax=ax[1], showfliers=showfliers,
                )
    ax[0].set(xlabel=r'sim. expected block length $\alpha$', )
    ax[1].set(xlabel=r'sim. expected block length $\alpha$',
              ylabel=r'deletion rate $\hat{\rho}_I$', title=r'Independent Deletion Model ($\alpha = 1$)')

    ax[1].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red')

    # ax[1].legend(title='Test result')
    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
                )
    # sns.stripplot(data=df, x='loss model param.', y='sim. loss rate', ax=ax[2], hue=None, marker='_', facecolor='red',
    #               edgecolor=None, size=100)
    for i, val in enumerate(unique_sim_lr):
        ax[2].axhline(val, i * 1 / len(unique_sim_lr), (i + 1) * 1 / len(unique_sim_lr), color='red')
    ax[2].set(xlabel=r'sim. expected block length $\alpha$',
              ylabel=r'deletion rate $\hat{\rho}_B$')  # , title=r'Block loss model: initialization rate'
    # ax[2].legend(labels=[r'simulation $\rho_B$'], )
    # ax[2].legend(title='Test result')
    sns.boxplot(data=df_plot, x='loss model param.', y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers, )

    unique_sim_alpha = pd.unique(df['loss model param.'].dropna())
    for i, val in enumerate(unique_sim_alpha):
        ax[3].axhline(val, i * 1 / len(unique_sim_alpha), (i + 1) * 1 / len(unique_sim_alpha), color='red')
    # unique_sim_alpha = pd.unique(df['loss model param.'])
    # for i, val in enumerate(unique_sim_alpha):
    #     ax[3].axhline(val, i*1/len(unique_sim_alpha), (i + 1)*1/len(unique_sim_alpha), color='red')
    ax[3].set(xlabel=r'sim. expected block length $\alpha$',
              ylabel=r'expected block length $\hat{\alpha}$')  # , title=r'Block loss model: expected block length'
    # ax[3].legend(title='Test result')
    plt.figtext(0.76, 0.95, 'Block Deletion Model', ha='center', va='center')

    df_plot.loc[:, 'loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    sns.boxplot(data=df_plot, x='loss model param.', y='loss_rate_1_times_alpha_1',
                hue=None, ax=ax[4], showfliers=showfliers, )
    ax[4].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red')
    ax[4].set(xlabel=r'sim. expected block length $\alpha$',
              ylabel=r'deletion rate per spacer $\hat{\rho}_B \cdot \hat{\alpha}$')

    df_plot.loc[:, 'gain_rate_1'] = df_plot['loss_rate_1_times_alpha_1'] * df_plot['avg array length']
    unique_sim_gain_rate = df_plot['sim. gain rate'].median()
    sns.boxplot(data=df_plot, x='loss model param.', y='gain_rate_1', hue=None, ax=ax[5], showfliers=showfliers, )
    ax[5].axhline(unique_sim_gain_rate, 0, len(unique_sim_alpha), color='red')
    ax[5].set(xlabel=r'sim. expected block length $\alpha$', ylabel=r'gain rate $\hat{\theta}_B$')
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_estimated_parameters_vs_x(df, save_path, showfliers=True, x='sim. loss rate'):
    """
    :param x:
    :param df:
    :param save_path:
    :param showfliers:
    :return:
    """

    unique_sim_lr = pd.unique(df['sim. loss rate'])
    fig, ax = plt.subplots(1, 6, figsize=(35, 5))
    sns.countplot(data=df, x=x, hue='test result', ax=ax[0], )
    # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    df_plot = df[df['all max loss lengths'].astype(bool)]
    ax[0].set(title='Test results')
    ax[0].legend(title='Test result', labels=['IDM', 'BDM'])
    sns.boxplot(data=df_plot, x=x, y='loss_rate_0', hue=None, ax=ax[1], showfliers=showfliers,
                )
    ax[0].set(xlabel=x, )
    ax[1].set(xlabel=x,
              ylabel=r'deletion rate $\hat{\rho}_I$', title=r'Independent Deletion Model ($\alpha = 1$)')

    # ax[1].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red')

    # ax[1].legend(title='Test result')
    sns.boxplot(data=df_plot, x=x, y='loss_rate_1', hue=None, ax=ax[2], showfliers=showfliers,
                )
    # sns.stripplot(data=df, x='loss model param.', y='sim. loss rate', ax=ax[2], hue=None, marker='_', facecolor='red',
    #               edgecolor=None, size=100)
    # for i, val in enumerate(unique_sim_lr):
    #     ax[2].axhline(val, i * 1 / len(unique_sim_lr), (i + 1) * 1 / len(unique_sim_lr), color='red')
    ax[2].set(xlabel=x,
              ylabel=r'deletion rate $\hat{\rho}_B$')  # , title=r'Block loss model: initialization rate'
    # ax[2].legend(labels=[r'simulation $\rho_B$'], )
    # ax[2].legend(title='Test result')
    sns.boxplot(data=df, x=x, y='alpha_1', hue=None, ax=ax[3], showfliers=showfliers, )

    # unique_sim_alpha = pd.unique(df['loss model param.'].dropna())
    # for i, val in enumerate(unique_sim_alpha):
    #     ax[3].axhline(val, i * 1 / len(unique_sim_alpha), (i + 1) * 1 / len(unique_sim_alpha), color='red')
    # unique_sim_alpha = pd.unique(df['loss model param.'])
    # for i, val in enumerate(unique_sim_alpha):
    #     ax[3].axhline(val, i*1/len(unique_sim_alpha), (i + 1)*1/len(unique_sim_alpha), color='red')
    ax[3].set(xlabel=x,
              ylabel=r'expected block length $\hat{\alpha}$')  # , title=r'Block loss model: expected block length'
    # ax[3].legend(title='Test result')
    plt.figtext(0.76, 0.95, 'Block Deletion Model', ha='center', va='center')

    df_plot.loc[:, 'loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    sns.boxplot(data=df_plot, x=x, y='loss_rate_1_times_alpha_1',
                hue=None, ax=ax[4], showfliers=showfliers, )
    # ax[4].axhline(max(unique_sim_lr), 0, len(unique_sim_lr), color='red')
    ax[4].set(xlabel=x,
              ylabel=r'deletion rate per spacer $\hat{\rho}_B \cdot \hat{\alpha}$')

    df_plot.loc[:, 'gain_rate_1'] = df_plot['loss_rate_1_times_alpha_1'] * df_plot['avg array length']
    unique_sim_gain_rate = df_plot['sim. gain rate'].median()
    sns.boxplot(data=df_plot, x=x, y='gain_rate_1', hue=None, ax=ax[5], showfliers=showfliers, )
    # ax[5].axhline(unique_sim_gain_rate, 0, len(unique_sim_alpha), color='red')
    ax[5].set(xlabel=x, ylabel=r'gain rate $\hat{\theta}_B$')
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_avg_array_length(df, save_path):
    sns.histplot(data=df, x='avg array length', stat='count', binwidth=1)
    plt.tight_layout(pad=3.0)
    plt.savefig(save_path)
    plt.show()


def plot_orientation_stats(df, save_path):
    df = df.explode('array orientation', ignore_index=True)
    fig, ax = plt.subplots(1, 4, figsize=(20, 5))
    sns.countplot(data=df, x='loss model param.', hue='recommend reversing array', ax=ax[0], )
    # ax[0].set_xticklabels(ax[0].get_xticks(), rotation=90)
    ax[0].set(title='Recommended reversion')
    # ax[0].legend(title='Recommended reversion', labels=['IDM', 'BDM'])
    sns.countplot(data=df, x='loss model param.', hue='array orientation', ax=ax[1], )
    ax[1].set(title='simulated orientation by reversing (2) or not (1)')
    mask = df['recommend reversing array'].values

    def f(x):
        if x == 1:
            return 2
        if x == 2:
            return 1

    df['predicted orientation'] = [f(o) if m else o for o, m in zip(df['array orientation'].values, mask)]
    sns.countplot(data=df, x='loss model param.', hue='predicted orientation', ax=ax[2], )
    ax[2].set(title='Predicted orientation by our method, (1) is correct')

    # df['orientation accuracy'] = [True if sim == pred else False
    #                               for sim, pred in
    #                               zip(df['array orientation'].values, df['predicted orientation'].values)]
    # sns.countplot(data=df, x='loss model param.', hue='orientation accuracy', ax=ax[3], )
    # ax[3].set(title='orientation accuracy')
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_orientation_likelihood_distribution(df, save_path, log=True):
    # lh_0_key = ''
    # reversed_lh_0_key = ''
    if log:
        lh_key = 'ln_lh_1'
        reversed_lh_key = 'reversed_ln_lh_1'
        lh_ratio = 'ln_lh_1 - reversed_ln_lh_1'
    else:
        lh_key = 'lh_1'
        reversed_lh_key = 'reversed_lh_1'
        lh_ratio = 'lh_1 / reversed_lh_1'
    fig, ax = plt.subplots(1, 3, figsize=(20, 5))
    sns.histplot(data=df, x=lh_key, ax=ax[0])
    sns.histplot(data=df, x=reversed_lh_key, ax=ax[1])
    sns.histplot(data=df, x=lh_ratio, hue='recommend reversing array', multiple='stack', ax=ax[2], bins=50)
    fig.suptitle('log likelihood distributions (orientation)' if log else 'likelihood distributions (orientation)')
    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def paper_sim_based_on_real_parameters(df, save_path, showfliers=False, stripplot_size=3, log_transform=False):
    df_plot = df.copy()

    df_plot = df_plot[df_plot['all max loss lengths'].astype(bool)]
    df_plot['loss_rate_1_times_alpha_1'] = df_plot['loss_rate_1'] * df_plot['alpha_1']
    df_plot['sim_loss_rate_times_alpha'] = df_plot['sim. loss rate'] * df_plot['loss model param.']
    df_plot['delta_loss_rate_1'] = df_plot['loss_rate_1'] - df_plot['sim. loss rate']
    df_plot['delta_alpha_1'] = df_plot['alpha_1'] - df_plot['loss model param.']
    df_plot['delta_loss_rate_1_times_alpha_1'] = (df_plot['loss_rate_1_times_alpha_1'] -
                                                  df_plot['sim_loss_rate_times_alpha'])

    if log_transform:
        df_plot['delta_loss_rate_1_times_alpha_1'] = np.log10(df_plot['delta_loss_rate_1_times_alpha_1'])
        df_plot['delta_loss_rate_1'] = np.log10(df_plot['delta_loss_rate_1'])
        df_plot['delta_alpha_1'] = np.log10(df_plot['delta_alpha_1'])

    fig = plt.figure(figsize=(15, 10))
    gs = GridSpec(1, 3, figure=fig)
    # title = fig.add_subplot(gs[0, :])
    ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])]

    sns.violinplot(data=df_plot,
                   x='sim loss model', y='delta_loss_rate_1_times_alpha_1',
                   ax=ax[0], showfliers=showfliers, inner='quartile', )
    # ylims = ax[0].get_ylim()
    sns.stripplot(data=df_plot, x='sim loss model', y='delta_loss_rate_1_times_alpha_1',
                  ax=ax[0], jitter=True, size=stripplot_size, dodge=True, legend=None,
                  edgecolor='gray', linewidth=1)
    ax[0].set(xlabel=r'$\rho_B \cdot \alpha$', ylabel=r'$\Delta \rho_{B} \cdot \alpha$')

    sns.violinplot(data=df_plot,
                   x='sim loss model', y='delta_loss_rate_1',
                   ax=ax[1], showfliers=showfliers, inner='quartile', )
    # ylims = ax[0].get_ylim()
    sns.stripplot(data=df_plot, x='sim loss model', y='delta_loss_rate_1',
                  ax=ax[1], jitter=True, size=stripplot_size, dodge=True, legend=None,
                  edgecolor='gray', linewidth=1)
    ax[1].set(xlabel=r'$\rho_B$', ylabel=r'$\Delta \rho_{B}$')

    sns.violinplot(data=df_plot,
                   x='sim loss model', y='delta_alpha_1',
                   ax=ax[2], showfliers=showfliers, inner='quartile', )
    # ylims = ax[0].get_ylim()
    sns.stripplot(data=df_plot, x='sim loss model', y='delta_alpha_1',
                  ax=ax[2], jitter=True, size=stripplot_size, dodge=True, legend=None,
                  edgecolor='gray', linewidth=1)
    ax[2].set(xlabel=r'$\alpha$', ylabel=r'$\Delta \alpha$')

    fig.tight_layout(pad=1.0)
    plt.show()
    fig.savefig(save_path)


def plot_fes_les_loss_freq(df, save_path, max_considered_distance=10, use_global=False):
    df_plot = df.copy()
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
    fig.suptitle('Global FES/LES/LS deletion frequencies'
                 if use_global else 'FES/LES/LS deletion frequencies')
    ax[0].bar(x=ls_fes_labels, height=ls_fes_freq, width=0.8, align='center')
    ax[0].set_xlabel('Positions to the right of FES')
    ax[0].set_ylabel('Deletion frequencies')
    ax[0].set_title('FES+ deletion frequencies')

    ax[1].bar(x=ls_les_labels[::-1], height=ls_les_freq[::-1], width=0.8, align='center')
    ax[1].set_xlabel('Positions to the left of LES')
    ax[1].set_ylabel('Deletion frequencies')
    ax[1].set_title('LES- deletion frequencies')

    if not use_global:
        ax[2].bar(x=ls_ls_labels[::-1], height=ls_ls_freq[::-1], width=0.8, align='center')
        ax[2].set_xlabel('Positions to the left of LS')
        ax[2].set_ylabel('Deletion frequencies')
        ax[2].set_title('LS- deletion frequencies')

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


def plot_paper_fs_ls_loss_freq(df, save_path, max_considered_distance=10, binning=False, sim=False, interlaced=False):
    df_plot = df.copy()

    sns.set_context("paper", font_scale=2)

    if binning:
        bins = [0, 15, 25, np.infty]
        names = ['<15', '15 - 25', '>25']
        df_plot['avg array length bins'] = pd.cut(df_plot['avg array length'], bins, labels=names)
        df_plot = df_plot[df_plot['avg array length bins'] == '15 - 25']

    fs_key = 'sim fs +m presence' if sim else 'fs +m presence'
    ls_key = 'sim ls -m presence' if sim else 'ls -m presence'

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
        ax.set_title('simulations')
    else:
        palette = sns.color_palette().as_hex()
        palette = [palette[1] if (i == 0) else palette[0] for i in range(max_considered_distance + 1)]
        # palette = ['tab:orange' if (i == 0) else 'tab:blue' for i in range(max_considered_distance + 1)]

        fig, ax = plt.subplots(1, 2, figsize=(20, 5))
        # fig.suptitle('simulation FS/LS deletion frequencies' if sim else 'FS/LS deletion frequencies')
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


if __name__ == '__main__':
    sim_on_real_parameters = False
    orientation = False
    sns.set_theme()
    # data_path = os.path.join('0_result_folder', 'sim_as_rec_new_old_lh_fct_exp')
    # data_path = os.path.join('0_result_folder', 'sim_determine_orientation')
    data_path = os.path.join('0_result_folder', '0_paper')
    # file_names = file_names + ['20_blm_1_5_len_15_leafs_16_gr_6_old_lh_fct_sim_as_rec_wo_unob']
    # file_names = ['22_blm_1_5_len_100_leafs_16_gr_80_new_lh_fct_sim_as_rec_w_unob']
    # file_names = file_names + ['23_blm_1_5_len_100_leafs_16_gr_160_new_lh_fct_sim_as_rec_w_unob']
    # file_names = file_names + ['24_blm_1_5_len_50_leafs_16_gr_80_new_lh_fct_sim_as_rec_w_unob']
    # file_names = ['8_frag_1_5_len_1000_leafs_16_gr_100']
    # file_names = ['9_test_flm']
    # file_names = ['11_based_on_real_data_block']
    # file_names = ['12_flm_gr_10_lr_001_tr_4-16']
    # file_names = ['13_flm_gr_10_lr_001_tr_16']
    # file_names = ['14_flm_gr_1000_lr_001_tr_16']
    # file_names = ['15_block_1_5_len_15_leafs_16_gr_6']
    # file_names = ['testing']
    # file_names = ['2_simulation_LRT_evaluation']
    # file_names = ['2_simulation_LRT_evaluation_nlh_sim_as_rec']
    # file_names = ['2_simulation_LRT_evaluation_rho_bias_correction']
    # file_names = ['2_simulation_LRT_evaluation_rho+alpha_bias_correction']
    # file_names = ['00_paper_simulation_optimal_estimation']
    file_names = ['00_paper_simulation_optimal_estimation_mean_length']
    # file_names = ['00_paper_real_para_based_sim_roc_curve']
    # file_names = ['00_paper_real_para_based_sim_roc_curve_wo_bias_old']
    # file_names = ['00_paper_sim_gdd', '00_paper_sim_gdd_bdy_corr']
    # file_names = ['00_paper_sim_sim_as_rec_wo_unobs', '00_paper_sim_sim_as_rec_w_unobs']
    # file_names = ['00_paper_simulation_optimal_estimation_small_trees']
    # file_names = ['00_paper_gdd_median_real_para_based_sim', '00_paper_gdd_sim_est_sim_est']
    # file_names = ['00_paper_gdd_mean_real_para_based_sim']
    # file_names = ['00_paper_gdd_median_real_para_based_sim_w_unob', '00_paper_gdd_median_real_para_based_sim_wo_unob',
    #               ]
    # file_names = ['00_paper_gdd_median_real_para_based_sim']
    # file_names = ['00_paper_gdd_mean_real_para_based_sim', ]
    # file_names = ['00_paper_simulation_optimal_estimation_high_spacer_turnover',]
    # file_names = ['00_paper_simulation_optimal_estimation_even_higher_spacer_turnover']
    # file_names = ['00_paper_sim_no_bias_corr']
    # file_names = ['00_paper_sim_precise_lh_fct']
    # file_names = ['00_paper_real_para_based_sim']
    # file_names = ['00_paper_simulation_optimal_estimation_bdy_corr']
    # file_names = ['00_paper_simulation_optimal_estimation_large_arrays']
    # file_names = ['2_simulation_LRT_evaluation_rho_bias_correction_only_idm']
    # file_names = ['2_simulation_LRT_evaluation_nlh_sim_as_rec_wo_unobserved']
    # file_names = ['2_simulation_LRT_evaluation_bidirectional']
    # file_names = ['2_simulation_LRT_evaluation_nlh_w_estimated_trees']
    # file_names = ['2_simulation_LRT_evaluation_nlh_w_estimated_trees_3_type_thetarho']
    # file_names = ['2_simulation_LRT_evaluation_bidirectional']
    # file_names = ['2_simulation_LRT_evaluation_nlh_sim_as_rec']
    protocol_name = '0_sim_rec_protocol.pkl'
    for name in file_names:
        # path = os.path.join('pictures', 'simulations', name, '0_exp_protocol.csv')
        # path = os.path.join('pictures', 'simulations', name, '0_exp_protocol.pkl')
        # path = os.path.join('pictures', 'simulations', name, '0_exp_protocol_wo_notconverged.pkl')
        path = os.path.join(data_path, name, protocol_name)
        save_path = os.path.join('1_result_graphs', '0_paper', name)
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        showfliers = False

        pd.set_option('display.max_columns', None)
        pd.set_option('display.max_rows', None)
        # df = pd.read_csv(path, index_col=0, converters={'relative loss positions': ast.literal_eval,
        #                                                 'sim. relative loss positions': ast.literal_eval,
        #                                                 'all max loss lengths': ast.literal_eval,
        #                                                 # 'not seen spacers': ast.literal_eval,
        #                                                 'all max loss lengths (normalized)': ast.literal_eval})
        df = pd.read_pickle(path)
        # print(df[max(df['alpha_1']) == df['alpha_1']])
        sns.set_style("whitegrid")
        sns.countplot(data=df, x='test result', )
        plt.tight_layout(pad=1.0)
        plt.savefig(os.path.join(save_path, '_'.join([name, 'test_results_count.pdf'])))
        plt.show()

        df_k = df[df['all max loss lengths'].astype(bool)]
        print('rho estimation', df_k['loss_rate_1'].mean())
        print('alpha estimation', df_k['alpha_1'].mean())
        print('rho median', df_k['loss_rate_1'].median())
        print('alpha median', df_k['alpha_1'].median())

        if orientation:
            ls_thresholds = list(range(0, 100, 1))
            plot_orientation_threshold_roc_curve(df, ls_thresholds,
                                                 os.path.join(save_path, '_'.join([name, 'roc_curve.pdf'])))
            plot_orientation_likelihood_distribution(df, os.path.join(save_path, '_'.join([name,
                                                                                           'orientation_ln_lh_dist.pdf'])),
                                                     log=True)

        if sim_on_real_parameters:
            paper_sim_based_on_real_parameters(df,
                                               os.path.join(save_path, '_'.join([name,
                                                                                 'paper_sim_based_on_real_parameters.pdf'])),
                                               showfliers=showfliers, )
        else:
            paper_plot_estimated_parameters(df,
                                            os.path.join(save_path, '_'.join([name, 'paper_estimated_parameters.pdf'])),
                                            showfliers=showfliers)

        plot_paper_fs_ls_loss_freq(df,
                                   os.path.join(save_path, '_'.join([name, 'paper_fs_ls_loss_freq_interlaced.pdf'])),
                                   max_considered_distance=5, binning=False, interlaced=True)
        plot_paper_fs_ls_loss_freq(df, os.path.join(save_path, '_'.join([name, 'paper_fs_ls_loss_freq.pdf'])),
                                   max_considered_distance=5, binning=False)
        plot_paper_fs_ls_loss_freq(df, os.path.join(save_path, '_'.join([name, 'paper_fs_ls_loss_freq_bin.pdf'])),
                                   max_considered_distance=5, binning=True)
        plot_paper_fs_ls_loss_freq(df, os.path.join(save_path, '_'.join([name, 'paper_fs_ls_loss_freq_sim.pdf'])),
                                   max_considered_distance=5, binning=False, sim=True)
        plot_paper_fs_ls_loss_freq(df, os.path.join(save_path, '_'.join([name, 'paper_fs_ls_loss_freq_bin_sim.pdf'])),
                                   max_considered_distance=5, binning=True, sim=True)

        plot_rel_pos_freq(df, os.path.join(save_path, '_'.join([name, 'rel_loss_frequencies.pdf'])), sim=False, )
        plot_rel_pos_freq(df, os.path.join(save_path, '_'.join([name, 'sim_rel_loss_frequencies.pdf'])), sim=True, )
        if 'sim. branch based relative loss positions' in df.columns:
            plot_rel_pos_freq(df, os.path.join(save_path, '_'.join([name, 'branch_based_sim_rel_loss_freq.pdf'])),
                              sim='branch_based_sim')
        # plot_orientation_likelihood_distribution(df, os.path.join(save_path, '_'.join([name,
        #                                                                                'orientation_lh_dist.pdf'])),
        #                                          log=False)

        # plot_rel_pos(df, os.path.join(save_path, '_'.join([name, 'rel_pos_hist.pdf'])))
        # plot_sim_rel_pos(df, os.path.join(save_path, '_'.join([name, 'sim_rel_pos_hist.pdf'])))
        plot_max_loss_lengths(df, os.path.join(save_path, '_'.join([name, 'max_loss_length.pdf'])))
        plot_max_loss_lengths(df, os.path.join(save_path, '_'.join([name, 'max_loss_length(100).pdf'])),
                              select_frag_par=1.0)
        # normalized to array length
        plot_max_loss_lengths_norm(df, os.path.join(save_path, '_'.join([name, 'max_loss_length_norm.pdf'])))
        plot_max_loss_lengths_norm(df, os.path.join(save_path, '_'.join([name, 'max_loss_length_norm(100).pdf'])),
                                   select_frag_par=1.0)
        plot_test_results_vs_frag_para(df, os.path.join(save_path, '_'.join([name, 'frag_loss_vs_test_result.pdf'])))
        plot_nb_gained_spacers(df, os.path.join(save_path, '_'.join([name, 'nb_gained_spacers.pdf'])))
        # plot_frag_para_bin_plot(df, os.path.join(save_path, '_'.join([name, 'frag_para_bin.pdf'])))
        plot_estimated_parameters(df, os.path.join(save_path, '_'.join([name, 'estimated_para.pdf'])),
                                  showfliers=showfliers)
        plot_avg_array_length(df, os.path.join(save_path, '_'.join([name, 'avg_array_len.pdf'])))
        plot_alternatives_estimated_parameters(df,
                                               os.path.join(save_path, '_'.join([name, 'diff_estimations_para.pdf'])),
                                               showfliers=showfliers)
        plot_estimated_parameters_vs_x(df, os.path.join(save_path, '_'.join([name, 'est_para_vs_sim_lr.pdf'])),
                                       x='sim. loss rate', showfliers=showfliers)
        plot_estimated_parameters_vs_x(df, os.path.join(save_path, '_'.join([name, 'est_para_vs_sim_gr.pdf'])),
                                       x='sim. gain rate', showfliers=showfliers)
        if orientation:
            plot_orientation_stats(df,
                                   os.path.join(save_path, '_'.join([name, 'orientation_stats.pdf'])))
