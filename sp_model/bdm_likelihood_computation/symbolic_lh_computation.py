import sympy
from sympy import Function, Derivative, exp, Eq, dsolve, simplify, init_printing, pprint, refine, Symbol, zeros, \
    summation, lambdify, symbols
import pickle
import dill
import os


def pickle_lambdified_lh_fct(save_path, lambdified_lh_fct):
    dill.settings['recurse'] = True
    with open(save_path + '.pickle', 'wb') as f:
        dill.dump(lambdified_lh_fct, f)


def load_pickled_lambdified_lh_fct(path):
    # Is this needed?
    # dill.settings['recurse'] = True
    with open(path, 'rb') as f:
        fct = dill.load(f)
    return fct


def load_lambdify_lh_fct(load_path, save_path=None, save_lambdified_lh_fct=True):
    """
    Lambdify a block deletion likelihood function. (This is for the pure death deletion model).
    :return:
    """
    lh_fct = load_p_from_pickle(load_path)
    a = symbols('a', positive=True)
    r, t = symbols('r t', nonnegative=True)
    lambdified_lh_fct = []
    # lambdifyed_ilm_lh_fct = []
    for lh in lh_fct[:, 0]:
        lambdified_lh_fct.append(lambdify((t, r, a), lh, 'numpy'))
        # sub_lh = lh.subs(a, 1)
        # lambdifyed_ilm_lh_fct.append(lambdify((t, r), sub_lh, 'numpy'))
    if save_lambdified_lh_fct:
        pickle_lambdified_lh_fct(save_path, lambdified_lh_fct)
    return lambdified_lh_fct


def solve_kolmogorov_backward_eq_birth_process(len_del, start_point):
    """
    param len_del: len_del is the deletion length + 1, i.e. the number of repeats.
    :param start_point:
    :return:
    """
    t = Symbol('t', nonnegative=True)
    r = Symbol('r', nonnegative=True)
    a = Symbol('a', positive=True)
    g_a = []
    p = zeros(len_del)
    for i in range(p.shape[0]):
        p[i, i] = sympy.exp(-r * t) ** (i + 1)
        g_a.append((1 - 1 / a) ** i)
    for j in range(2, p.shape[0] + 1):
        for i in range(j - 1, 0, -1):
            f = Function('f')
            rs = i * r / a * sum([g_a[k] * (p[k + i, j - 1] - f(t)) for k in range(0, j - i)])
            print(rs)
            print(i, j)
            deq = Eq(Derivative(f(t), t), rs)
            sol = dsolve(deq, f(t), ics={f(0): 0})
            p[i - 1, j - 1] = sol.args[1]
            print(p)
    return p


def solve_kolmogorov_backward_eq_death_process(len_del, log=True, save_path=None, cont_p=None):
    t = Symbol('t', nonnegative=True)
    r = Symbol('r', nonnegative=True)
    a = Symbol('a', positive=True)

    p = zeros(len_del + 1, 1)
    p_final = zeros(len_del + 1, 1)
    if cont_p is not None:
        for i in range(cont_p.shape[0]):
            p[i, 0] = cont_p[i, 0]
        start_value = cont_p.shape[0]
        print('Start value', start_value)
        print('Start vector', cont_p)
    else:
        start_value = 1
    g_a = []
    p[0, 0] = 1
    for i in range(p.shape[0]):
        g_a.append((1 - 1 / a) ** i)
    for i in range(start_value, p.shape[0]):
        f = Function('f')
        rs = r / a * sum([g_a[i - k - 1] * (k + 1) * p[k, 0] for k in range(0, i)]) - i * r * f(t)
        print('Differential equation: ', rs)
        deq = Eq(Derivative(f(t), t), rs)
        sol = dsolve(deq, f(t), ics={f(0): 0})
        p[i, 0] = sol.args[1]
        print('Solution: ', p[i, 0])
        if i % 1 == 0 and save_path is not None:
            save_p_to_pickle(save_path, p[:i, 0])
    if log:
        for i in range(p.shape[0]):
            p_final[i, 0] = simplify(sympy.log(p[i, 0]))
    else:
        for i in range(p.shape[0]):
            p_final[i, 0] = p[i, 0]
    return p_final


def apply_log(p):
    p_new = zeros(p.shape[0], 1)
    for i in range(p.shape[0]):
        print('i: ', i)
        print('Previous: ', p[i, 0])
        p_new[i, 0] = simplify(sympy.log(p[i, 0]))
        print('After ', p_new[i, 0])
    return p_new


def save_p_to_pickle(save_path, p):
    with open(save_path + '.pickle', 'wb') as f:
        pickle.dump(p, f)
    return


def load_p_from_pickle(save_path):
    with open(save_path, 'rb') as f:
        p = pickle.load(f)
    return p


if __name__ == '__main__':
    load_lambdify_lh_fct(os.path.join('sympy_bdm_lh_fct', '230329_death_lh_up_to_68.pickle'),
                         save_path=os.path.join('sympy_bdm_lh_fct',
                                                '230329_death_lh_up_to_68_lambdifyed'),
                         save_lambdified_lh_fct=True)
    # solve_kolmogorov_backward_eq(4)
    #
    # init_printing()
    # t = Symbol('t', nonnegative=True)
    # r = Symbol('r', nonnegative=True)
    # a = Symbol('a', nonnegative=True)
    # f = Function('f')
    # eq = Eq(Derivative(f(t), t), r / a * ((exp(-r * t)) ** 2 - f(t)))
    #
    # solution = dsolve(eq, f(t))
    # sim_solution = simplify(solution)
    # pprint(solution)
    # pprint(sim_solution)

    # p = solve_kolmogorov_backward_eq_birth_process(3, 0)
    # p = solve_kolmogorov_backward_eq_death_process(3, log=False)
    # pprint(p)
    # print('bla')
    # p = load_p_from_pickle(os.path.join('../../data', '0_lh', 'death_lh_up_to_100.pickle'))
    # print(len(p))
    # print('bla')
    # p_ilm = []
    # for i in range(p.shape[0] - 20):
    #     p_ilm.append(p[i, 0].subs(a ** 1, 1.0))
    #     print(i)
    #     pprint(p_ilm[i])
    # if i < 15:
    #     pprint(p_ilm[i])
    # pprint(p_ilm)
    # p = solve_kolmogorov_backward_eq_death_process(100, log=False,
    #                                                save_path=os.path.join('../../data', '0_lh',
    #                                                                       'death_lh_up_to_100_groot'),
    #                                                # cont_p=p
    #                                                cont_p=None
    #                                                )
    # p_log = apply_log(p)
    # save_p_to_pickle(os.path.join('additional_data', '0_lh', 'death_lh_up_to_100_groot'), p)

    # Depending on version of sympy/numpy/...? Might need to re-lambdify the function, or do it on loading instead ->
    # more expensive.
    load_lambdify_lh_fct(os.path.join('../../additional_data', '0_lh', '230329_death_lh_up_to_100_groot.pickle'),
                         save_path=os.path.join('../../additional_data', '0_lh',
                                                '230329_death_lh_up_to_100_groot_lambdifyed'),
                         save_lambdified_lh_fct=True)

    # p = load_p_from_pickle(os.path.join('data', '0_lh', 'testing_lh.pickle'))
    # a, r, t = symbols('a r t', nonnegative=True)
    # f = lambdify((t, r, a), p[10, 0])
    # print(p)
    # pprint(f)
    # pprint(p[10, 0])
    # pprint(p[0,0])
    # pprint(p[1,0])
    # pprint(p[2,0])
    # pprint(f(1, 0.1, 1))
