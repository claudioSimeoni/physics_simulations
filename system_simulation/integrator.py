import numpy as np


# with Hamiltonian Mechanics the state is assumed to be q = state[:n], p = state[n:]
class DynamicSystem:
    def __init__(self, time, state, number_of_bodies, f = None, hp = None, hq = None): # TODO: handling errors if I call a not defined function 
        self.time = time
        self.state = state
        self.n = number_of_bodies
        self.f = f
        self.hp = hp
        self.hq = hq


class Integrator:
    def __init__(self, h):
        self.h = h

    def update(self, system):
        raise NotImplementedError
    

class RK4(Integrator):
    def update(self, system):
        dt, t, y, f = self.h, system.time, system.state, system.f
        k1 = f(t, y)
        k2 = f(t + dt/2, y + dt/2 * k1)
        k3 = f(t + dt/2, y + dt/2 * k2)
        k4 = f(t + dt, y + dt * k3)
        system.t, system.state = t + dt, y + dt/6 * (k1 + 2 * k2 + 2 * k3 + k4)
        

class Leapfrog(Integrator):
    def update(self, system):
        q, p, hp, hq, dt = system.state[:system.n], system.state[system.n:], system.hp, system.hq, self.h
        pm = p - dt / 2 * hq(q)
        q = q + dt * hp(pm)
        p = pm - dt / 2 * hq(q)

# TODO: implementing other integrators

