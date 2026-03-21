# with Hamiltonian the state is assumed to be q = state[:npc], p = state[npc:]
class DynamicSystem:
    def __init__(self, state, state_vector_dimension, space_dimension):
        self.time = 0
        self.state = state
        self.dim = state_vector_dimension
        self.npc = state_vector_dimension // 2 # number of positional coordinates
        self.nb = state_vector_dimension // (2 * space_dimension) # number of bodies


class Integrator:
    def __init__(self, h):
        self.h = h

    def update(self, system):
        raise NotImplementedError
    

class RK4(Integrator):
    def update(self, system, f):
        dt, t, y, f = self.h, system.time, system.state, f
        k1 = f(t, y)
        k2 = f(t + dt/2, y + dt/2 * k1)
        k3 = f(t + dt/2, y + dt/2 * k2)
        k4 = f(t + dt, y + dt * k3)
        system.time, system.state = t + dt, y + dt/6 * (k1 + 2 * k2 + 2 * k3 + k4)


class ExplicitEuler(Integrator):
    def update(self, system, Hp, Hq):
        q, p, dt = system.state[:system.npc], system.state[system.npc:], self.h
        q[:], p[:] = q + dt * Hp(p), p + dt * -Hq(q) # contextual assignment otherwise we would need another variable
        system.time += dt
        

class SymplecticEuler(Integrator):
    def update(self, system, Hp, Hq):
        q, p, dt = system.state[:system.npc], system.state[system.npc:], self.h
        p[:] = p + dt * -Hq(q)
        q[:] = q + dt * Hp(p)
        system.time += dt
    

class VerletStormer(Integrator):
    def update(self, system, Hp, Hq):
        q, p, dt = system.state[:system.npc], system.state[system.npc:], self.h
        pm = p - dt / 2 * Hq(q)
        q[:] = q + dt * Hp(pm)
        p[:] = pm - dt / 2 * Hq(q)
        system.time += dt

# TODO: implementing other integrators

