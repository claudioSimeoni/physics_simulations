import numpy as np

n = 3
G = 1

planets_mass = 1
dq0 = np.array([0, (-3 ** 0.25) / 2, (3 ** 0.25) / 2, (3 ** -0.25), (-3 ** -0.25) / 2, (-3 ** -0.25) / 2])

m = np.ones(n) * planets_mass
mm = m[:, None] * m[None, :]
q = np.array([1, -1/2, -1/2, 0, (3 ** 0.5) / 2, (-3 ** 0.5) / 2])
p = np.concatenate((dq0[:n] / m, dq0[n:] / m))

sizes = 100 * np.arctan(m / 10)
soft = 0.01 * np.mean(np.abs(q)) # constant for softening

ANIMATION_LENGTH = 10
FRAMES_PER_SECOND = 60
FRAMES = ANIMATION_LENGTH * FRAMES_PER_SECOND
ANIMATION_INTERVAL = 1 / FRAMES_PER_SECOND
UPDATES_PER_FRAME = 1
DT = 1 / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)