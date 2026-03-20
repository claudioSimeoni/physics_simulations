import numpy as np
from astroquery.jplhorizons import Horizons
from astropy.time import Time
import astropy.units as u


planet_ids = ['0', '1', '2', '3', '4', '5', '6', '7', '8']
names = ["Sun", "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"]


epoch = Time.now().jd

R = []
V = []

print("Fetching data from NASA JPL Horizons...")

for p_id in planet_ids:
    # 'location=@ssb' sets the origin to the Solar System Barycenter
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


e = 1
G = 6.67430e-20  # km^3 kg^-1 s^-2
GM_VALUES = [
    132712440041.93938, # Sun
    22031.78,           # Mercury
    324858.59,          # Venus
    398600.435,         # Earth
    42828.37,           # Mars
    126686531.91,       # Jupiter
    37931184.3,         # Saturn
    5793939,            # Uranus
    6835100             # Neptune
]
PLANETS_MASSES = [gm / G for gm in GM_VALUES]


# ------------------------------------------------------------------------------
#                              exported values
# ------------------------------------------------------------------------------


n = 9
m = np.array(PLANETS_MASSES)
mm = m[:, None] * m[None, :]
q = np.array([R[i][0] for i in range(n)] +
             [R[i][1] for i in range(n)] + 
             [R[i][2] for i in range(n)])
p = np.concatenate((np.array([V[i][0] for i in range(n)]) * m * e, 
                    np.array([V[i][1] for i in range(n)]) * m * e, 
                    np.array([V[i][2] for i in range(n)]) * m * e))


# non physical constants
color = [
    "yellow",       # Sun
    "gray",         # Mercury
    "orange",       # Venus
    "blue",         # Earth
    "red",          # Mars
    "orange",       # Jupiter
    "lightblue",    # Saturn
    "cyan",         # Uranus
    "blueviolet"    # Neptune
]
sizes = [
    696340/2000,    # Sun
    2440/2000,      # Mercury
    6052/2000,      # Venus
    6371/2000,      # Earth
    3390/2000,      # Mars
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
UPDATES_PER_FRAME = 1
DAYS_PER_SECOND = 356
TIME_PASSED_IN_A_SECOND = 60 * 60 * 24 * DAYS_PER_SECOND
DT = TIME_PASSED_IN_A_SECOND / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)