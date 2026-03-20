from .data.solar_system_NASA import *
from collections import deque
from .. import integrator
from .. import animator

# this is project specific, contains
#                                   1. the equations
#                                   2. the update function of the gui
#                                   3. all the parameters
# 

class Simulator:
    def __init__(self, n, m, mm, soft, sizes, color, integ, system):
        self.n = n
        self.m = m
        self.mm = mm
        self.soft = soft
        self.sizes = sizes
        self.color = color
        self.integ = integ
        self.system = system

    def H(self): 
        q, p = self.system.state[:3 * n], self.system.state[3 * n:]
        T1 = np.sum((p[:n] ** 2) / (2 * m))
        T2 = np.sum((p[n:2 * n] ** 2) / (2 * m))
        T3 = np.sum((p[2 * n:] ** 2) / (2 * m))

        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:2 * n, None] - q[None, n:2 * n]
        z_diff = q[2 * n:, None] - q[None, 2 * n:]
        module = (x_diff ** 2 + y_diff ** 2 + z_diff ** 2) ** 0.5
        np.fill_diagonal(module, np.inf)

        U = np.sum(- mm / module) * G / 2

        return T1 + T2 + T3 + U

    def Hp(self, p):
        return np.concatenate((p[:n] / self.m, p[n:2 * n] / self.m, p[2 * n:] / self.m))

    def Hq(self, q):
        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:2 * n, None] - q[None, n:2 * n]
        z_diff = q[2 * n:, None] - q[None, 2 * n:]
        module = (x_diff ** 2 + y_diff ** 2 + z_diff ** 2 + self.soft ** 2) ** 1.5
        np.fill_diagonal(module, np.inf)
        return np.concatenate((np.sum((self.mm * G * x_diff / module), axis=1),
                               np.sum((self.mm * G * y_diff / module), axis=1),
                               np.sum((self.mm * G * z_diff / module), axis=1)))

    def step(self):
        self.integ.update(self.system, self.Hp, self.Hq)


class BodiesPlotting:
    def __init__(self, xlim, ylim, zlim, sim, proj):
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim
        self.nb = len(sim.system.state) // 6 # number of bodies
        self.coord = sim.system.state[:3 * self.nb]
        self.sizes = sim.sizes
        self.color = sim.color
        self.proj = proj

    def setup(self, ax):
        self.ax = ax
        ax.set(xlim3d=(-self.xlim, self.xlim), xlabel='X')
        ax.set(ylim3d=(-self.ylim, self.ylim), ylabel='Y')
        ax.set(zlim3d=(-self.zlim, self.zlim), zlabel='Z')
        self.ax.set_facecolor("white")
        self.ax.set_aspect("equal")
        self.ax.grid()

        nb = self.nb
        self.graphics = ax.scatter(self.coord[:nb], self.coord[nb:2 * nb], self.coord[2 * nb:], 
                                   s=self.sizes, c=self.color)

    def update(self):
        xs = self.coord[:self.nb]
        ys = self.coord[self.nb:2 * self.nb]
        zs = self.coord[2 * self.nb:3 * self.nb]

        self.graphics._offsets3d = (xs, ys, zs)
        return self.graphics,


# class TracePlotting:
#     def __init__(self, xlim, ylim, zlim, tracelength, sim, proj):
#         self.xlim = xlim
#         self.ylim = ylim
#         self.zlim = zlim
#         self.nb = len(sim.system.state) // 4
#         self.coord = sim.system.state[:2 * self.nb]
#         self.sizes = sim.sizes
#         self.tracelength = tracelength
#         self.proj = proj

#         self.pos = [[deque(maxlen=tracelength) for _ in range(self.nb)], 
#                     [deque(maxlen=tracelength) for _ in range(self.nb)]]

#     def setup(self, ax):
#         self.ax = ax
#         self.ax.set_xlim(-self.xlim, self.xlim)
#         self.ax.set_ylim(-self.ylim, self.ylim)
#         self.ax.set_facecolor("white")
#         self.ax.set_aspect("equal")
#         self.ax.grid()

#         self.graphics = ax.scatter(self.coord[:self.nb], self.coord[self.nb:], 
#                                    s=self.sizes, c="black")

#         from matplotlib.lines import Line2D
#         self.trace_lines = []
#         for _ in range(self.nb):
#             line = Line2D([], [], color="black", lw=2, alpha=0.7)
#             self.ax.add_line(line)
#             self.trace_lines.append(line)

#     def update(self):
#         for i in range(self.nb):
#             self.pos[0][i].append(self.coord[i])
#             self.pos[1][i].append(self.coord[i + self.nb])

#         for i in range(self.nb):
#             self.trace_lines[i].set_data(self.pos[0][i], self.pos[1][i])

#         self.graphics.set_offsets(np.column_stack((self.coord[:self.nb], 
#                                                    self.coord[self.nb:])))
#         return [self.graphics] + self.trace_lines


class EnergyPlotting:
    def __init__(self, sim, proj):
        self.t = []
        self.energy = []
        self.sim = sim
        self.proj = proj

    def setup(self, ax):
        self.ax = ax
        self.ax.grid()

        self.line, = self.ax.plot(self.t, self.energy, label="Energy")

    def update(self):
        self.t.append(self.sim.system.time)
        self.energy.append(self.sim.H())
        self.line.set_data(self.t, self.energy)

        self.ax.relim()
        self.ax.autoscale_view()
        return self.line,
    

if __name__ == "__main__":
    integ = integrator.VerletStormer(DT)
    system = integrator.DynamicSystem(np.concatenate((q, p)), 6 * n, 3)
    sim = Simulator(n, m, mm, soft, sizes, color, integ, system)

    lim = np.max(np.abs(q)) * 3

    bp = BodiesPlotting(lim, lim, lim, sim, "3d")
    # tp = TracePlotting(lim, lim, lim, 100, sim, "3d")
    ep = EnergyPlotting(sim, None)
    ani = animator.Animator(ANIMATION_LENGTH, FRAMES_PER_SECOND, UPDATES_PER_FRAME,
                            sim, [bp])
    
    ani.run()