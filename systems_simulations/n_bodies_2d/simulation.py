from .configuration.infinite_symbol import *
from collections import deque
from .. import integrator
from .. import animator


'''
The Simulator stores everything regarding the physical simulation: the DynamicSystem, the update function 
associated with the specific DynamicSystem and the Integrator, that is used in the step function to update 
the state of the DynamicSystem. 

The Animator from animator.py then uses the Simulator to update the DynamicSystem, and the Plotting (e.g.
BodiesPlotting) to update the matplotlib gui. In this way to add a new plot associated with a Simulator,
you can simply define a new Plotting object. (TODO: It might make sense to make a Plotting from which the others inherit?)
'''


class Simulator:
    def __init__(self, n, m, mm, soft, integ, system):
        self.n = n
        self.m = m
        self.mm = mm
        self.soft = soft
        self.integ = integ
        self.system = system

    def H(self): 
        n, m, mm = self.n, self.m, self.mm
        q, p = self.system.state[:2 * n], self.system.state[2 * n:]
        T1 = np.sum((p[:n] ** 2) / (2 * m))
        T2 = np.sum((p[n:] ** 2) / (2 * m))

        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:, None] - q[None, n:]
        module = (x_diff ** 2 + y_diff ** 2 + self.soft ** 2) ** 0.5
        np.fill_diagonal(module, np.inf)
        U = np.sum(-mm / module) * G / 2

        return T1 + T2 + U

    def Hp(self, p):
        n, m = self.n, self.m
        return np.concatenate((p[:n] / m, p[n:] / m))

    def Hq(self, q):
        n, mm = self.n, self.mm
        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:, None] - q[None, n:]
        module = (x_diff ** 2 + y_diff ** 2 + self.soft ** 2) ** 1.5
        np.fill_diagonal(module, np.inf)
        return np.concatenate((np.sum((mm * G * x_diff / module), axis=1),
                               np.sum((mm * G * y_diff / module), axis=1)))

    def step(self):
        self.integ.update(self.system, self.Hp, self.Hq)


class BodiesPlotting:
    def __init__(self, xlim, ylim, sizes, color, sim, proj):
        self.xlim = xlim
        self.ylim = ylim
        self.nb = len(sim.system.state) // 4 # number of bodies
        self.coord = sim.system.state[:2 * self.nb]
        self.sizes = sizes
        self.color = color
        self.proj = proj

    def setup(self, ax):
        self.ax = ax
        self.ax.set_xlim(-self.xlim, self.xlim)
        self.ax.set_ylim(-self.ylim, self.ylim)
        self.ax.set_facecolor("white")
        self.ax.set_aspect("equal")
        self.ax.grid()

        nb = self.nb
        self.graphics = ax.scatter(self.coord[:nb], self.coord[nb:], 
                                   s=self.sizes, c=self.color)

    def update(self):
        self.graphics.set_offsets(np.column_stack((self.coord[:self.nb], self.coord[self.nb:])))
        return self.graphics,


class TracePlotting:
    def __init__(self, xlim, ylim, sizes, color, tracelength, sim, proj):
        self.xlim = xlim
        self.ylim = ylim
        self.nb = len(sim.system.state) // 4
        self.coord = sim.system.state[:2 * self.nb]
        self.sizes = sizes
        self.color = color
        self.tracelength = tracelength
        self.proj = proj

        # each deque stores [x, y]
        self.traces = [deque(maxlen=tracelength) for _ in range(self.nb)]

    def setup(self, ax):
        self.ax = ax
        self.ax.set_xlim(-self.xlim, self.xlim)
        self.ax.set_ylim(-self.ylim, self.ylim)
        self.ax.set_facecolor("white")
        self.ax.set_aspect("equal")
        self.ax.grid()

        self.graphics = ax.scatter(self.coord[:self.nb], self.coord[self.nb:], s=self.sizes, c=self.color)

        # list of ax elements, each of them is a trace (a line plot)
        self.trace_lines = []
        for _ in range(self.nb):
            line, = ax.plot([], [], color="black", lw=2, alpha=0.1)
            self.trace_lines.append(line)

    def update(self):
        x = self.coord[:self.nb]
        y = self.coord[self.nb:]

        for i, trace in enumerate(self.traces):
            trace.append([x[i], y[i]])

        for line, trace in zip(self.trace_lines, self.traces):
            pts = np.array(trace)
            line.set_data(pts[:, 0], pts[:, 1])

        self.graphics.set_offsets(np.column_stack((x, y)))
        return [self.graphics] + self.trace_lines


class EnergyPlotting:
    def __init__(self, sim, proj):
        self.t = []
        self.energy = []
        self.sim = sim
        self.proj = proj

    def setup(self, ax):
        self.ax = ax
        self.ax.grid()
        self.line, = self.ax.plot(self.t, self.energy, label="Energy", color='red', lw=2)
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Energy (J)")
        self.ax.set_title("Energy")
        self.ax.legend()

    def update(self):
        self.t.append(self.sim.system.time)
        self.energy.append(self.sim.H())
        self.line.set_data(self.t, self.energy)

        self.ax.relim()
        self.ax.autoscale_view()
        return self.line,
    

if __name__ == "__main__":
    integ = integrator.VerletStormer(DT)
    system = integrator.DynamicSystem(np.concatenate((q, p)), 4 * n, 2)
    sim = Simulator(n, m, mm, soft, integ, system)

    bp = BodiesPlotting(np.max(np.abs(q)) * 3, np.max(np.abs(q)) * 3, sizes, color, sim, None)
    tp = TracePlotting(np.max(np.abs(q)) * 3, np.max(np.abs(q)) * 3, sizes, color, 300, sim, None)
    ep = EnergyPlotting(sim, None)
    ani = animator.Animator("white", ANIMATION_LENGTH, FRAMES_PER_SECOND, UPDATES_PER_FRAME,
                            sim, [tp, ep]) # in the list insert any Plotting objects you want to plot
    
    ani.run()