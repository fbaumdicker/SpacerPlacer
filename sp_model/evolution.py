import numpy as np
import scipy


class GTR:
    def __init__(self, gain_rate, loss_rate, k=1):
        self.gain_rate = gain_rate
        self.loss_rate = loss_rate
        self.k = k

        self.q = self.set_q()
        self.stationary_p = self.set_stationary_prob()

    def modify_rates(self, new_gain_rate, new_loss_rate, new_k=1):
        self.gain_rate = new_gain_rate
        self.loss_rate = new_loss_rate
        self.k = new_k

        self.q = self.set_q()
        self.stationary_p = self.set_stationary_prob()
        return

    def set_q(self):
        q = np.zeros((2, 2))
        q[0, 0] = - self.gain_rate / self.k
        q[0, 1] = self.gain_rate / self.k
        q[1, 0] = self.loss_rate
        q[1, 1] = -self.loss_rate
        return q

    def set_stationary_prob(self):
        norm = self.gain_rate / self.k + self.loss_rate
        ratio_0 = self.loss_rate / norm
        ratio_1 = self.gain_rate / self.k / norm
        prob = np.zeros((2, 2))
        prob[:, 0] = ratio_0
        prob[:, 1] = ratio_1
        return prob

    def prob(self, t):
        qt = scipy.linalg.expm(self.q * t)
        return np.maximum(0, qt)

    def evolve(self, profile, t, return_log=False):
        qt = self.prob(t)
        res = profile.dot(qt)
        return np.log(res) if return_log else res

    def prob_loss(self, t):
        return max(0, 1 - np.exp(-self.loss_rate * t))

    def propagate_profile(self, profile, t, return_log=False):
        qt = self.prob(t).T
        res = profile.dot(qt)
        return np.log(res) if return_log else res

    def stationary_prof(self, profile, return_log=False):
        prob = self.stationary_p
        res = profile.dot(prob)
        return np.log(res) if return_log else res
