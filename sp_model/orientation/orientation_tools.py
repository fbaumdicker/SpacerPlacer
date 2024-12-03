import pandas as pd


def compare_likelihoods_for_orientation(protocol, protocol_reversed_raw, protocol_boring=None,
                                        protocol_boring_reversed=None,
                                        decision_boundary=0, dict_trees=None):
    """

    :param dict_trees:
    :param protocol:
    :param protocol_reversed_raw:
    :param protocol_boring:
    :param protocol_boring_reversed:
    :param decision_boundary:
    :return:
    """
    if dict_trees:
        dict_tree_forward = dict_trees['forward']
        dict_tree_reversed = dict_trees['reversed']
    else:
        dict_tree_forward = None
        dict_tree_reversed = None

    if protocol_boring is not None:
        if not protocol_boring.empty:
            protocol = pd.concat([protocol, protocol_boring])
        if not protocol_boring_reversed.empty:
            protocol_reversed_raw = pd.concat([protocol_reversed_raw, protocol_boring_reversed])

    protocol_reversed = protocol_reversed_raw.rename(columns={'lh_idm': 'reversed_lh_idm',
                                                              'ln_lh_idm': 'reversed_ln_lh_idm',
                                                              'lh_bdm': 'reversed_lh_bdm',
                                                              'ln_lh_bdm': 'reversed_ln_lh_bdm'})

    new_protocol = pd.merge(protocol, protocol_reversed[['reversed_lh_idm',
                                                         'reversed_ln_lh_idm',
                                                         'reversed_lh_bdm',
                                                         'reversed_ln_lh_bdm']],
                            left_index=True, right_index=True, how='left')

    new_protocol['lh_idm / reversed_lh_idm'] = new_protocol['lh_idm'] / new_protocol['reversed_lh_idm']
    new_protocol['ln_lh_idm - reversed_ln_lh_idm'] = new_protocol['ln_lh_idm'] - new_protocol['reversed_ln_lh_idm']
    new_protocol['lh_bdm / reversed_lh_bdm'] = new_protocol['lh_bdm'] / new_protocol['reversed_lh_bdm']
    new_protocol['ln_lh_bdm - reversed_ln_lh_bdm'] = new_protocol['ln_lh_bdm'] - new_protocol['reversed_ln_lh_bdm']

    new_protocol['recommend reversing array'] = [decision_fct(v, decision_boundary=decision_boundary)
                                                 for v in new_protocol['ln_lh_bdm - reversed_ln_lh_bdm']]

    # final_columns = list(new_protocol.columns) + ['our predicted orientation']
    oriented_protocol = protocol.copy()
    protocol_boring_index = set() if protocol_boring is None or protocol_boring.empty else set(protocol_boring.index)
    protocol_boring_reversed_index = set() if protocol_boring_reversed is None or protocol_boring_reversed.empty \
        else set(protocol_boring_reversed.index)

    ls_boring = []
    ls_recommendation = []
    ls_pred_orient = []
    ls_reversed_lh_0, ls_reversed_lh_1, ls_reversed_ln_lh_0, ls_reversed_ln_lh_1 = [], [], [], []
    dict_trees_final = {}
    for i, n in enumerate(protocol.index):
        if new_protocol.loc[n, 'recommend reversing array'] == True:
            protocol.loc[n, :] = protocol_reversed_raw.loc[n, :]
            ls_reversed_lh_0.append(new_protocol.loc[n, 'lh_idm'])
            ls_reversed_lh_1.append(new_protocol.loc[n, 'lh_bdm'])
            ls_reversed_ln_lh_0.append(new_protocol.loc[n, 'ln_lh_idm'])
            ls_reversed_ln_lh_1.append(new_protocol.loc[n, 'ln_lh_bdm'])
            ls_pred_orient.append('Reverse')
            ls_recommendation.append(True)
            if n in protocol_boring_reversed_index:
                ls_boring.append(n)

            if dict_trees:
                dict_trees_final[n] = dict_tree_reversed[n]
        else:
            if n in protocol_boring_index:
                ls_boring.append(n)
            ls_reversed_lh_0.append(protocol_reversed_raw.loc[n, 'lh_idm'])
            ls_reversed_lh_1.append(protocol_reversed_raw.loc[n, 'lh_bdm'])
            ls_reversed_ln_lh_0.append(protocol_reversed_raw.loc[n, 'ln_lh_idm'])
            ls_reversed_ln_lh_1.append(protocol_reversed_raw.loc[n, 'ln_lh_bdm'])
            ls_pred_orient.append('Forward' if not new_protocol.loc[n, 'recommend reversing array'] else 'ND')
            ls_recommendation.append(new_protocol.loc[n, 'recommend reversing array'])

            if dict_trees:
                dict_trees_final[n] = dict_tree_forward[n]
    oriented_protocol['recommend reversing array'] = ls_recommendation
    oriented_protocol['predicted orientation'] = ls_pred_orient
    oriented_protocol['reversed_lh_idm'] = ls_reversed_lh_0
    oriented_protocol['reversed_lh_bdm'] = ls_reversed_lh_1
    oriented_protocol['reversed_ln_lh_idm'] = ls_reversed_ln_lh_0
    oriented_protocol['reversed_ln_lh_bdm'] = ls_reversed_ln_lh_1
    oriented_protocol['lh_idm / reversed_lh_idm'] = new_protocol['lh_idm / reversed_lh_idm']
    oriented_protocol['ln_lh_idm - reversed_ln_lh_idm'] = new_protocol['ln_lh_idm - reversed_ln_lh_idm']
    oriented_protocol['lh_bdm / reversed_lh_bdm'] = new_protocol['lh_bdm / reversed_lh_bdm']
    oriented_protocol['ln_lh_bdm - reversed_ln_lh_bdm'] = new_protocol['ln_lh_bdm - reversed_ln_lh_bdm']
    oriented_protocol = oriented_protocol.drop(ls_boring, axis=0, errors='ignore')

    # print(oriented_protocol)
    return new_protocol, oriented_protocol, dict_trees_final


def decision_fct(ln_lh_diff, decision_boundary=0):
    decision = 'nd'
    if ln_lh_diff < -decision_boundary:
        decision = True
    elif ln_lh_diff > decision_boundary:
        decision = False
    return decision
