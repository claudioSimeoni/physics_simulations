from data.random_100_bodies import *
import integrator

# this is project specific, contains
#                                   1. the equations
#                                   2. the update function of the gui
#                                   3. all the parameters

class Simulator:
    def __init__(self, n, m, mm, soft):
        self.n = n
        self.m = m
        self.mm = mm
        self.soft = soft

    def hp(p):
        return np.concatenate((p[:n] / m, p[n:] / m))

    def hq(q):
        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:, None] - q[None, n:]
        module = (x_diff ** 2 + y_diff ** 2 + soft ** 2) ** 1.5
        np.fill_diagonal(module, np.inf)
        return np.concatenate((np.sum((mm * x_diff / module), axis=1), np.sum((mm * y_diff / module), axis=1)))
    




