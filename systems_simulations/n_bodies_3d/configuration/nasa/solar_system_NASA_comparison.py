import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.time import Time
import astropy.units as u


planet_ids = ['10', '1', '2', '399', '301', '4', '5', '6', '7', '8']
names = ["Sun", "Mercury", "Venus", "Earth", "Moon", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]

n = len(planet_ids)
G = 6.67430e-20  # km^3 kg^-1 s^-2

GM_VALUES = [
    132712440041.93938, # Sun
    22031.78,           # Mercury
    324858.59,          # Venus
    398600.4354,        # Earth
    4902.800,           # Moon
    42828.37,           # Mars
    126686531.91,       # Jupiter
    37931184.3,         # Saturn
    5793939,            # Uranus
    6835100             # Neptune
]

PLANETS_MASSES = [gm / G for gm in GM_VALUES]

R = []
V = []


def retrieve_solar_system_data(date):
    global R, V
    R = []
    V = []

    epoch = Time(date, scale='utc').tdb.jd
    
    for p_id in planet_ids:
        obj = Horizons(id=p_id, location='@ssb', epochs=epoch)
        vec = obj.vectors()
        
        R.append([
            vec['x'].quantity.to(u.km).value[0], 
            vec['y'].quantity.to(u.km).value[0], 
            vec['z'].quantity.to(u.km).value[0]
        ])
        
        V.append([
            vec['vx'].quantity.to(u.km/u.s).value[0], 
            vec['vy'].quantity.to(u.km/u.s).value[0], 
            vec['vz'].quantity.to(u.km/u.s).value[0]
        ])

    # print("\n", date, ": ")
    # for i in range(0, n):
    #     print(names[i], ": \n", "R = ", R[i], "\n", "V = ", V[i])


m = np.array(PLANETS_MASSES)
mm = m[:, None] * m[None, :]
q = np.array([])
p = np.array([])


def compute_initial_positions():
    global q, p
    q = np.array([R[i][0] for i in range(n)] +
                 [R[i][1] for i in range(n)] + 
                 [R[i][2] for i in range(n)])
    
    p = np.concatenate((np.array([V[i][0] for i in range(n)]) * m, 
                        np.array([V[i][1] for i in range(n)]) * m, 
                        np.array([V[i][2] for i in range(n)]) * m))


def revert_pq(q, p):
    R = [[q[i], q[i+n], q[i+2*n]] for i in range(n)]
    return R


# animation constants
DAY_FRACTION_FOR_SIMULATION = 1
DT = 60 * 60 * 24 / DAY_FRACTION_FOR_SIMULATION