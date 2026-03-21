import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.time import Time
import astropy.units as u


planet_ids = ['10', '5', '6', '7', '8']
names = ["Sun", "Jupiter", "Saturn", "Uranus", "Neptune"]

epoch = Time("2026-03-20", scale='utc').tdb.jd

R = []
V = []

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


e = 0
G = 6.67430e-20  # km^3 kg^-1 s^-2
GM_VALUES = [
    132712440041.93938, # Sun
    126686531.91,       # Jupiter
    37931184.3,         # Saturn
    5793939,            # Uranus
    6835100             # Neptune
]
PLANETS_MASSES = [gm / G for gm in GM_VALUES]


# ------------------------------------------------------------------------------
#                              exported values
# ------------------------------------------------------------------------------


n = 5
m = np.array(PLANETS_MASSES)
mm = m[:, None] * m[None, :]
q = np.array([R[i][0] for i in range(n)] +
             [R[i][1] for i in range(n)] + 
             [R[i][2] for i in range(n)])
p = np.concatenate((np.array([V[i][0] for i in range(n)]) * m * (1 + e), 
                    np.array([V[i][1] for i in range(n)]) * m * (1 + e), 
                    np.array([V[i][2] for i in range(n)]) * m * (1 + e)))


# non physical constants
color = [
    "yellow",       # Sun
    "orange",       # Jupiter
    "lightblue",    # Saturn
    "cyan",         # Uranus
    "blueviolet"    # Neptune
]
sizes = [
    696340/2000,    # Sun
    69911/2000,     # Jupiter
    58232/2000,     # Saturn
    25362/2000,     # Uranus
    24622/2000      # Neptune
]
soft = 0.01 * np.mean(np.abs(q))


# animation constants
ANIMATION_LENGTH = 10
FRAMES_PER_SECOND = 60
FRAMES = ANIMATION_LENGTH * FRAMES_PER_SECOND
ANIMATION_INTERVAL = 1 / FRAMES_PER_SECOND
UPDATES_PER_FRAME = 6
DAYS_PER_SECOND = 356
TIME_PASSED_IN_A_SECOND = 60 * 60 * 24 * DAYS_PER_SECOND
DT = TIME_PASSED_IN_A_SECOND / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)