import numpy as np
from scipy.spatial.transform import Rotation
import random


n = 100
G = 1
PLANETS_MASSES = [100000] + [1] * (n - 1)
STAR_MASS = PLANETS_MASSES[0]
e = 0 # eccentricity of orbits, can even be a np.array()


def rotate_orbits(n, R, V):
    for i in range(n):
        axis = np.array([random.random() for _ in range(3)])
        axis = axis / np.linalg.norm(axis)
        angle = random.random() * 2 * np.pi

        rot = Rotation.from_rotvec(angle * axis)
        R[i] = rot.apply(np.array(R[i])).tolist()
        V[i] = rot.apply(np.array(V[i])).tolist()


# 1. R, V initialization through LIST OF INITIAL PERIODS / MODULES OF SPEED AND COORDINATES / CAN ROTATE VECTORS WITH THE ROTATE

MAX_T = 5
T = np.array([0] + [MAX_T * random.random() for _ in range(n - 1)])
MOD_R = (((T ** 2) * STAR_MASS * G) / (4 * (np.pi ** 2))) ** (1 / 3)
MOD_V = np.concatenate((np.zeros((1, )), (STAR_MASS * G / MOD_R[1:]) ** (1 / 2)))

R = [[MOD_R[i].item(), 0, 0] for i in range(n)]
V = [[0, MOD_V[i].item(), 0] for i in range(n)]

rotate_orbits(n, R, V)


# 2. R, V initialization through LIST OF COORDINATES / SPEEDS

# R = [[1, 2, 3], 
#      [4, 5, 6],
#      [7, 8, 9]]

# V = [[1, 2, 3], 
#      [4, 5, 6],
#      [7, 8, 9]]


# ------------------------------------------------------------------------------
#                              exported values
# ------------------------------------------------------------------------------


n
m = np.array(PLANETS_MASSES)
mm = m[:, None] * m[None, :]
q = np.array([R[i][0] for i in range(n)] +
             [R[i][1] for i in range(n)] + 
             [R[i][2] for i in range(n)])
p = np.concatenate((np.array([V[i][0] for i in range(n)]) * m * (1 + e), 
                    np.array([V[i][1] for i in range(n)]) * m * (1 + e), 
                    np.array([V[i][2] for i in range(n)]) * m * (1 + e)))


# non physical constants
color = [np.random.rand(3) for _ in range(n)] 
# color = "black"
sizes = 100 * np.arctan(m / 10)
soft = 0.01 * np.mean(np.abs(q))


# animation constants
ANIMATION_LENGTH = 10
FRAMES_PER_SECOND = 60
FRAMES = ANIMATION_LENGTH * FRAMES_PER_SECOND
ANIMATION_INTERVAL = 1 / FRAMES_PER_SECOND
UPDATES_PER_FRAME = 1
TIME_PASSED_IN_A_SECOND = 1
DT = TIME_PASSED_IN_A_SECOND / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)