import scipy.stats
import numpy as np


def test_significance_ratio_chi2(value_to_test, significance_value, n=1):
    quantile = scipy.stats.chi2.isf(significance_value, n)
    return True if value_to_test > quantile else False, quantile


def test_ls_significance_ratio_chi2(ls_values_to_test, significance_value, n=1):
    ls_quantiles = []
    ls_test_results = []
    for val in ls_values_to_test:
        r, q = test_significance_ratio_chi2(val, significance_value, n=n)
        ls_quantiles.append(q)
        ls_test_results.append(r)
    return ls_test_results, ls_quantiles


def geom(alpha, k):
    if alpha == 1:
        if k == 1:
            return 1
        else:
            return 0
    p = 1 / alpha
    return (1 - p) ** (k - 1) * p


def np_rate_exponential(rate, *args, **kwargs):
    if rate == 0:
        return np.random.exponential(np.infty, *args, **kwargs)
    return np.random.exponential(1 / rate, *args, **kwargs)


def p_stationary_frag_length(n, rate):
    prod = 1
    p = np.ones(n + 1)
    for i in range(n + 1):
        val = (i + 1) * (i + 2)
        prod *= val / (2 * rate) + 1
        p[i] = val / (2 * rate * prod)
        # print(val, prod, p[i])
    return p


def sample_stationary_frag_length(sample_size, rate, dist_limit=1000):
    """
    Samples the stationary distribution of the fragment model. This distribution is created empirically, the dist_limit
    determines the precision (in theory dist_limit = \infty). (see Kupczok et al.)
    :param sample_size:
    :param rate:
    :param dist_limit:
    :return:
    """
    p = p_stationary_frag_length(dist_limit, rate)
    return np.random.choice(range(dist_limit + 1), sample_size, p=p / sum(p), replace=True)


def ev_stationary_frag_length(rate, dist_limit=1000):
    p = p_stationary_frag_length(dist_limit, rate)
    return sum([n * p[n] for n in range(len(p))])


def var_stationary_frag_length(rate, dist_limit=1000):
    p = p_stationary_frag_length(dist_limit, rate)
    ev = sum([n * p[n] for n in range(len(p))])
    ev_2 = sum([n ** 2 * p[n] for n in range(len(p))])
    return ev_2 - ev ** 2


def binomial_thinning_alpha(alpha, p_observed):
    return p_observed * (alpha - 1) + 1


def bionmial_thinning_inverse_alpha(alpha_observed, p_observed):
    return (alpha_observed - 1) / p_observed + 1





import matplotlib.pyplot as plt
import seaborn as sns

# xs = list(range(20))
# f = p_stationary_frag_length(len(xs)-1, 1)
# sample = sample_stationary_frag_length(100000, 10, dist_limit=1000000)
# fig = plt.figure()
# sns.histplot(sample, stat='probability')
# # plt.plot(xs, f)
# plt.show()
# print(exp_stationary_frag_length(6/0.2))
# print(exp_stationary_frag_length(6/0.3))
