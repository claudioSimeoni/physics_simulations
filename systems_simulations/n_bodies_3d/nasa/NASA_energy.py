from ..configuration.nasa import solar_system_NASA_comparison as nasa
from ... import integrator
from datetime import datetime, timedelta, date
import numpy as np
import matplotlib.pyplot as plt


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
        

    def H(self): 
        n, m, mm = self.n, self.m, self.mm
        q, p = self.system.state[:3 * n], self.system.state[3 * n:]
        T1 = np.sum((p[:n] ** 2) / (2 * m))
        T2 = np.sum((p[n:2 * n] ** 2) / (2 * m))
        T3 = np.sum((p[2 * n:] ** 2) / (2 * m))

        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:2 * n, None] - q[None, n:2 * n]
        z_diff = q[2 * n:, None] - q[None, 2 * n:]
        module = (x_diff ** 2 + y_diff ** 2 + z_diff ** 2) ** 0.5
        np.fill_diagonal(module, np.inf)

        U = np.sum(- mm / module) * nasa.G / 2

        return T1 + T2 + T3 + U

    def Hp(self, p):
        n, m = self.n, self.m
        return np.concatenate((p[:n] / m, p[n:2 * n] / m, p[2 * n:] / m))

    def Hq(self, q):
        n, mm = self.n, self.mm
        x_diff = q[:n, None] - q[None, :n]
        y_diff = q[n:2 * n, None] - q[None, n:2 * n]
        z_diff = q[2 * n:, None] - q[None, 2 * n:]
        module = (x_diff ** 2 + y_diff ** 2 + z_diff ** 2) ** 1.5
        np.fill_diagonal(module, np.inf)
        return np.concatenate((np.sum((mm * nasa.G * x_diff / module), axis=1),
                               np.sum((mm * nasa.G * y_diff / module), axis=1),
                               np.sum((mm * nasa.G * z_diff / module), axis=1)))

    def step(self):
        self.integ.update(self.system, self.Hp, self.Hq)


if __name__ == "__main__":
    days_of_simulation = 365 * 20
    points_to_plot = 10
    days_between_plotting = days_of_simulation // points_to_plot
    planet = 1

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

    e01 = sim1.H()
    e02 = sim2.H()
    e03 = sim3.H()
    e04 = sim4.H()

    energy1 = [0]
    energy2 = [0]
    energy3 = [0]
    energy4 = [0]

    days = [0]

    for i in range(points_to_plot):
        days.append((i + 1) * days_between_plotting)
        print((i + 1) / points_to_plot * 100, "%")
        for j in range(days_between_plotting * nasa.DAY_FRACTION_FOR_SIMULATION):
            sim1.step()
            sim2.step()
            sim3.step()
            sim4.step()
        
        energy1.append(np.abs((sim1.H() - e01) / e01) * 100)
        energy2.append(np.abs((sim2.H() - e02) / e02) * 100)
        energy3.append(np.abs((sim3.H() - e03) / e03) * 100)
        energy4.append(np.abs((sim4.H() - e04) / e04) * 100)

    print()
    print("Data with timestep ", 1 / nasa.DAY_FRACTION_FOR_SIMULATION, ":\n")
    print("Explicit Euler's error: ", max(energy1), "%")
    print("Symplectic Euler's error: ", max(energy2), "%")
    print("Verlet Stormer's error: ", max(energy3), "%")
    print("RK4's error: ", max(energy4), "%")

    ax.plot(days, energy1, label="Explicit Euler")
    ax.plot(days, energy2, label="Symplectic Euler")
    ax.plot(days, energy3, label="Verlet Stormer")
    ax.plot(days, energy4, label="RK4")

    ax.text(0.80, 0.95, f"Timestep = {1 / nasa.DAY_FRACTION_FOR_SIMULATION} days", transform=ax.transAxes)

    ax.set_ylabel(f"Energy deviation [%]")
    ax.set_xlabel("Time [days]")

    ax.legend()
    plt.show()

