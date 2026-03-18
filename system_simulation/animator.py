import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation


'''
- TODO: possibilities on what to plot: 
    1. Trajectories (with camera on center of mass)
    2. Shadows
    3. Energy
    4. Angular Momentum
    ...

- how will this work:
    the Animator takes the Elements to plot (object that describe each ax), when updating it 
    makes the system step on, and after that it plots the graphics with the updates of each 
    Element
'''


class Animator:
    def __init__(self, length, fps, updates_per_frame, simulator, visual_elements):
        self.length = length
        self.fps = fps
        self.updates_per_frame = updates_per_frame
        self.sim = simulator
        self.elems = visual_elements

        n = len(self.elems)
        self.fig = plt.figure()
        self.axs = []

        for i, el in enumerate(self.elems):
            ax = self.fig.add_subplot(1, n, i + 1, projection=el.proj) # TODO: maybe make more flexible the positioning
            el.setup(ax)
            self.axs.append(ax)

    def update(self, frame):
        for _ in range(self.updates_per_frame):
            self.sim.step()

        graphics = []
        for elem in self.elems:
            graphics.extend(elem.update())

        return graphics

    def run(self):
        self.ani = animation.FuncAnimation(
            fig = self.fig, 
            func = self.update, 
            frames = self.fps * self.length, 
            interval = 1000 / self.fps, # ms
            blit = False,
        )

        plt.show()