from .configuration.nasa import solar_system_NASA_comparison as nasa
from .. import integrator
from datetime import datetime, timedelta, date
import numpy as np


'''
This script (along with data/solar_system_NASA_comparison.py) is used only to retrieve NASA data of
`days_of_simulation` ago, computing today position according to this simulator, and comparing it with
today actual data from NASA. TODO: making the comparison easier to read.
'''


class Simulator:
    def __init__(self, n, m, mm, integ, system):
        self.n = n
        self.m = m
        self.mm = mm
        self.integ = integ
        self.system = system

    def Hp(self, p):
        return np.concatenate((
            p[:self.n] / self.m,
            p[self.n:2 * self.n] / self.m,
            p[2 * self.n:] / self.m
        ))

    def Hq(self, q):
        x_diff = q[:self.n, None] - q[None, :self.n]
        y_diff = q[self.n:2 * self.n, None] - q[None, self.n:2 * self.n]
        z_diff = q[2 * self.n:, None] - q[None, 2 * self.n:]
        module = (x_diff ** 2 + y_diff ** 2 + z_diff ** 2) ** 1.5
        np.fill_diagonal(module, np.inf)

        return np.concatenate((
            np.sum((self.mm * nasa.G * x_diff / module), axis=1),
            np.sum((self.mm * nasa.G * y_diff / module), axis=1),
            np.sum((self.mm * nasa.G * z_diff / module), axis=1)
        ))

    def step(self):
        self.integ.update(self.system, self.Hp, self.Hq)


if __name__ == "__main__":
    days_of_simulation = 365 * 10

    # subtracting days_of_simulation to today
    d = datetime.strptime(str(date.today()), "%Y-%m-%d").date()
    beginning_date = str(d - timedelta(days=days_of_simulation))

    # initial retrieving
    nasa.retrieve_solar_system_data(beginning_date)
    nasa.compute_initial_positions()

    # simulating
    integ = integrator.VerletStormer(nasa.DT)
    system = integrator.DynamicSystem(np.concatenate((nasa.q, nasa.p)), 6 * nasa.n, 3)
    sim = Simulator(nasa.n, nasa.m, nasa.mm, integ, system)

    for i in range(days_of_simulation * nasa.DAY_FRACTION_FOR_SIMULATION):
        sim.step()

    # reverting in R / V format and printing my data
    n = nasa.n
    nasa.q = sim.system.state[:3 * n]
    nasa.p = sim.system.state[3 * n:]
    R = [[nasa.q[i], nasa.q[i+n], nasa.q[i+2*n]] for i in range(n)]
    V = [[nasa.p[i] / (nasa.m[i]),
          nasa.p[i+n] / (nasa.m[i]),
          nasa.p[i+2*n] / (nasa.m[i])] for i in range(n)]
    
    print("\nmy data: ")

    for i in range(0, 9):
        print(nasa.names[i], " : \n", "R = ", R[i], "\n", "V = ", V[i])

    # retrieving today data for comparison
    nasa.retrieve_solar_system_data(str(date.today()))
    nasa.compute_initial_positions()
