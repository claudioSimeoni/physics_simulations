from ..configuration.nasa import solar_system_NASA_comparison as nasa
from ... import integrator
from datetime import datetime, timedelta, date
import numpy as np
import matplotlib.pyplot as plt
import sys


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
    days_of_simulation = 365 * 2
    points_to_plot = 10
    days_between_plotting = days_of_simulation // points_to_plot
    planet = 9

    # subtracting days_of_simulation to today
    d = datetime.strptime(str(date.today()), "%Y-%m-%d").date()
    beginning_date = str(d - timedelta(days=days_of_simulation))

    # initial retrieving
    nasa.retrieve_solar_system_data(beginning_date)
    nasa.compute_initial_positions()

    # simulating
    system1 = integrator.DynamicSystem(np.concatenate((nasa.q, nasa.p)), 6 * nasa.n, 3)
    system2 = integrator.DynamicSystem(np.concatenate((nasa.q, nasa.p)), 6 * nasa.n, 3)
    system3 = integrator.DynamicSystem(np.concatenate((nasa.q, nasa.p)), 6 * nasa.n, 3)
    system4 = integrator.DynamicSystem(np.concatenate((nasa.q, nasa.p)), 6 * nasa.n, 3)
    
    sim1 = Simulator(nasa.n, nasa.m, nasa.mm, integrator.ExplicitEuler(nasa.DT), system1)
    sim2 = Simulator(nasa.n, nasa.m, nasa.mm, integrator.SymplecticEuler(nasa.DT), system2)
    sim3 = Simulator(nasa.n, nasa.m, nasa.mm, integrator.VerletStormer(nasa.DT), system3)
    sim4 = Simulator(nasa.n, nasa.m, nasa.mm, integrator.RK4(nasa.DT), system4)

    # plotting
    fig, ax = plt.subplots()

    dist1 = [0]
    dist2 = [0]
    dist3 = [0]
    dist4 = [0]

    days = [0]
    for i in range(points_to_plot):
        days.append((i + 1) * days_between_plotting)
        print((i + 1) / points_to_plot * 100, "%", file=sys.stderr)
        for j in range(days_between_plotting * nasa.DAY_FRACTION_FOR_SIMULATION):
            sim1.step()
            sim2.step()
            sim3.step()
            sim4.step()
        
        d = datetime.strptime(beginning_date, "%Y-%m-%d").date()
        new_date = str(d + timedelta(days=(i + 1) * days_between_plotting))
        nasa.retrieve_solar_system_data(new_date)

        R1 = nasa.revert_pq(sim1.system.state[:3 * nasa.n], sim1.system.state[3 * nasa.n:])
        R2 = nasa.revert_pq(sim2.system.state[:3 * nasa.n], sim2.system.state[3 * nasa.n:])
        R3 = nasa.revert_pq(sim3.system.state[:3 * nasa.n], sim3.system.state[3 * nasa.n:])
        R4 = nasa.revert_pq(sim4.system.state[:3 * nasa.n], sim4.system.state[3 * nasa.n:])

        d1 = np.sqrt((R1[planet][0] - nasa.R[planet][0]) ** 2 +
                     (R1[planet][1] - nasa.R[planet][1]) ** 2 +
                     (R1[planet][2] - nasa.R[planet][2]) ** 2)
        d2 = np.sqrt((R2[planet][0] - nasa.R[planet][0]) ** 2 +
                     (R2[planet][1] - nasa.R[planet][1]) ** 2 +
                     (R2[planet][2] - nasa.R[planet][2]) ** 2)
        d3 = np.sqrt((R3[planet][0] - nasa.R[planet][0]) ** 2 +
                     (R3[planet][1] - nasa.R[planet][1]) ** 2 +
                     (R3[planet][2] - nasa.R[planet][2]) ** 2)
        d4 = np.sqrt((R4[planet][0] - nasa.R[planet][0]) ** 2 +
                     (R4[planet][1] - nasa.R[planet][1]) ** 2 +
                     (R4[planet][2] - nasa.R[planet][2]) ** 2)
        
        dist1.append(d1)
        dist2.append(d2)
        dist3.append(d3)
        dist4.append(d4)

    print()
    print(nasa.names[planet], "'s data with timestep ", 1 / nasa.DAY_FRACTION_FOR_SIMULATION, ":\n")
    print("Explicit Euler's error: ", dist1[len(dist1) - 1], "km")
    print("Symplectic Euler's error: ", dist2[len(dist2) - 1], "km")
    print("Verlet Stormer's error: ", dist3[len(dist3) - 1], "km")
    print("RK4's error: ", dist4[len(dist4) - 1], "km\n")

    ax.plot(days, dist1, label="Explicit Euler")
    ax.plot(days, dist2, label="Symplectic Euler")
    ax.plot(days, dist3, label="Verlet Stormer")
    ax.plot(days, dist4, label="RK4")

    ax.text(0.80, 0.95, f"Timestep = {1 / nasa.DAY_FRACTION_FOR_SIMULATION} days", transform=ax.transAxes)

    ax.set_ylabel(f"Error [km]")
    ax.set_xlabel("Time [days]")

    ax.legend()
    plt.show()