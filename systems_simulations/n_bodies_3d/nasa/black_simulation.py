from ..configuration.nasa.full_solar_system_NASA import *
from collections import deque
from ... import integrator
from ... import animator
import time


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
        q, p = self.system.state[:3 * n], self.system.state[3 * n:]
        T1 = np.sum((p[:n] ** 2) / (2 * m))
        T2 = np.sum((p[n:2 * n] ** 2) / (2 * m))
        T3 = np.sum((p[2 * n:] ** 2) / (2 * m))

        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:2 * n, None] - q[None, n:2 * n]
        z_diff = q[2 * n:, None] - q[None, 2 * n:]
        module = (x_diff ** 2 + y_diff ** 2 + z_diff ** 2 + self.soft ** 2) ** 0.5
        np.fill_diagonal(module, np.inf)

        U = np.sum(- mm / module) * G / 2

        return T1 + T2 + T3 + U

    def Hp(self, p):
        n, m = self.n, self.m
        return np.concatenate((p[:n] / m, p[n:2 * n] / m, p[2 * n:] / m))

    def Hq(self, q):
        n, mm = self.n, self.mm
        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:2 * n, None] - q[None, n:2 * n]
        z_diff = q[2 * n:, None] - q[None, 2 * n:]
        module = (x_diff ** 2 + y_diff ** 2 + z_diff ** 2 + self.soft ** 2) ** 1.5
        np.fill_diagonal(module, np.inf)
        return np.concatenate((np.sum((mm * G * x_diff / module), axis=1),
                               np.sum((mm * G * y_diff / module), axis=1),
                               np.sum((mm * G * z_diff / module), axis=1)))

    def step(self):
        self.integ.update(self.system, self.Hp, self.Hq)


class BodiesPlotting:
    def __init__(self, xlim, ylim, zlim, sizes, color, proj, start_time, sim):
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim
        self.nb = len(sim.system.state) // 6 # number of bodies
        self.coord = sim.system.state[:3 * self.nb]
        self.sizes = sizes
        self.color = color
        self.proj = proj

        # fps section
        self.start_time = start_time
        self.current_second = 0
        self.current_second_frames = 0
        self.fps = 0

    # if you want a traditional plot (no full black bg): set_facecolor("white"), grid(True), 
    # comment the 3 set_pane_color, set the texts color="white", (even the Animator bg="white" in main)
    def setup(self, ax):
        self.ax = ax
        ax.set(xlim3d=(-self.xlim, self.xlim), xlabel='X')
        ax.set(ylim3d=(-self.ylim, self.ylim), ylabel='Y')
        ax.set(zlim3d=(-self.zlim, self.zlim), zlabel='Z')
        ax.set_title("Galaxy")
        ax.set_facecolor("black")
        ax.set_aspect("equal")
        ax.grid(False)

        # comment this for white bg along with set_facecolor("white")
        ax.xaxis.set_pane_color((0,0,0,1))
        ax.yaxis.set_pane_color((0,0,0,1))
        ax.zaxis.set_pane_color((0,0,0,1))
        
        # this makes the plotting bigger
        ax.set_position([0, 0, 1, 1])

        # text
        self.time_text = ax.text2D(0.05, 0.95, '', transform=ax.transAxes, fontsize=10, color="white")
        self.fps_text = ax.text2D(0.95, 0.95, '', transform=ax.transAxes, fontsize=10, color="white")


        xs, ys, zs = self.coord[:self.nb], self.coord[self.nb:2*self.nb], self.coord[2*self.nb:]
        self.graphics = ax.scatter(xs, ys, zs, s=self.sizes, c=self.color)
        
    def update_fps(self):
        current_time = time.time()
        self.time_text.set_text(f"Time = {current_time - self.start_time}")

        if current_time - self.current_second >= 1:
            self.current_second = int(current_time)
            self.fps = self.current_second_frames
            self.current_second_frames = 0

        self.current_second_frames += 1
        self.fps_text.set_text(f"FPS = {self.fps}")

    def update(self):
        xs = self.coord[:self.nb]
        ys = self.coord[self.nb:2 * self.nb]
        zs = self.coord[2 * self.nb:3 * self.nb]

        self.update_fps()
        
        self.graphics._offsets3d = (xs, ys, zs)
        return self.graphics,


# TODO: maybe make this inherit?
class TracePlotting:
    def __init__(self, xlim, ylim, zlim, sizes, color, tracelength, proj, start_time, sim):
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim
        self.nb = len(sim.system.state) // 6 # number of bodies
        self.coord = sim.system.state[:3 * self.nb]
        self.sizes = sizes
        self.color = color
        self.tracelength = tracelength
        self.proj = proj
        self.sim = sim

        # each deque stores [x, y, z]
        self.traces = [deque(maxlen=tracelength) for _ in range(self.nb)]

        # fps tracking
        self.start_time = start_time
        self.current_second = 0
        self.current_second_frames = 0
        self.fps = 0

    def setup(self, ax):
        self.ax = ax
        ax.set(xlim3d=(-self.xlim, self.xlim), xlabel='X')
        ax.set(ylim3d=(-self.ylim, self.ylim), ylabel='Y')
        ax.set(zlim3d=(-self.zlim, self.zlim), zlabel='Z')
        ax.set_title("Galaxy")
        ax.set_facecolor("black")
        ax.set_aspect("equal")
        ax.grid(False)

        # comment this for white bg along with set_facecolor("white")
        ax.xaxis.set_pane_color((0,0,0,1))
        ax.yaxis.set_pane_color((0,0,0,1))
        ax.zaxis.set_pane_color((0,0,0,1))

        # this makes the plotting bigger
        ax.set_position([0, 0, 1, 1])

        # text
        self.time_text = ax.text2D(0.05, 0.95, '', transform=ax.transAxes, fontsize=10, color="white")
        self.fps_text = ax.text2D(0.95, 0.95, '', transform=ax.transAxes, fontsize=10, color="white")

        # list of ax elements, each of them is a trace (a line plot)
        self.trace_lines = []
        for _ in range(self.nb):
            line, = ax.plot([], [], [], color="white", lw=1.5, alpha=0.1)
            self.trace_lines.append(line)


        xs, ys, zs = self.coord[:self.nb], self.coord[self.nb:2*self.nb], self.coord[2*self.nb:]
        self.graphics = ax.scatter(xs, ys, zs, s=self.sizes, c=self.color)

    def update_fps(self):
        current_time = time.time()
        self.time_text.set_text(f"Time = {current_time - self.start_time}s")

        if current_time - self.current_second >= 1:
            self.current_second = int(current_time)
            self.fps = self.current_second_frames
            self.current_second_frames = 0

        self.current_second_frames += 1
        self.fps_text.set_text(f"FPS = {self.fps}")

    def update(self):
        xs = self.coord[:self.nb]
        ys = self.coord[self.nb:2*self.nb]
        zs = self.coord[2*self.nb:]

        self.update_fps()

        for i, trace in enumerate(self.traces):
            trace.append([xs[i], ys[i], zs[i]])

        for line, trace in zip(self.trace_lines, self.traces):
            if trace:
                pts = np.array(trace)
                line.set_data(pts[:, 0], pts[:, 1])
                line.set_3d_properties(pts[:, 2])

        self.graphics._offsets3d = (xs, ys, zs)
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
        self.energy_text = self.ax.text(0.05, 0.95, '', transform=self.ax.transAxes, fontsize=10)

    def update(self):
        self.t.append(self.sim.system.time)
        self.energy.append(self.sim.H())
        self.line.set_data(self.t, self.energy)

        self.energy_text.set_text(f"Energy = {self.energy[-1]}")

        self.ax.relim()
        self.ax.autoscale_view()
        return self.line,
    

if __name__ == "__main__":
    integrators = [integrator.ExplicitEuler(h=DT), integrator.SymplecticEuler(h=DT), 
                   integrator.VerletStormer(h=DT), integrator.RK4(h=DT)]
    integ_index = int(input('''
Select the integrator:
                            
0 = Explicit Euler
1 = Symplectic Euler
2 = Verlet Stormer
3 = RK4
                            
'''))
    
    integ = integrators[integ_index]
    system = integrator.DynamicSystem(state=np.concatenate((q, p)), state_vector_dimension=6 * n, space_dimension=3)
    sim = Simulator(n, m, mm, soft, integ, system)

    lim = np.max(np.abs(q)) * 3

    bp = BodiesPlotting(lim, lim, lim, sizes, color, proj="3d", start_time=time.time(), sim=sim)
    tp = TracePlotting(lim, lim, lim, sizes, color, tracelength=50, proj="3d", start_time=time.time(), sim=sim)
    ep = EnergyPlotting(sim, proj=None) # note that to plot this you gotta set everything to white otherwise looks awful
    ani = animator.Animator("black", ANIMATION_LENGTH, FRAMES_PER_SECOND, UPDATES_PER_FRAME,
                            sim, visual_elements=[tp]) # in the list insert any Plotting objects you want to plot
    
    ani.run()