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
e = 0
PLANETS_MASSES = [gm / G for gm in GM_VALUES]

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


# ------------------------------------------------------------------------------
#                              exported values
# ------------------------------------------------------------------------------


n = 10
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
    "gray",         # Mercury
    "orange",       # Venus
    "blue",         # Earth
    "silver",       # Moon
    "red",          # Mars
    "peru",         # Jupiter
    "khaki",        # Saturn
    "cyan",         # Uranus
    "royalblue"     # Neptune
]
sizes_scale = 1000
sizes = [
    696340/sizes_scale,    # Sun
    2440/sizes_scale,      # Mercury
    6052/sizes_scale,      # Venus
    6371/sizes_scale,      # Earth
    1737/sizes_scale,      # Moon
    3390/sizes_scale,      # Mars
    69911/sizes_scale,     # Jupiter
    58232/sizes_scale,     # Saturn
    25362/sizes_scale,     # Uranus
    24622/sizes_scale     # Neptune
]
soft = 0


# animation constants
ANIMATION_LENGTH = 10
FRAMES_PER_SECOND = 60
FRAMES = ANIMATION_LENGTH * FRAMES_PER_SECOND
ANIMATION_INTERVAL = 1 / FRAMES_PER_SECOND
UPDATES_PER_FRAME = 1
DAYS_PER_SECOND = 180
TIME_PASSED_IN_A_SECOND = 60 * 60 * 24 * DAYS_PER_SECOND
DT = TIME_PASSED_IN_A_SECOND / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)