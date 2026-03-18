import numpy as np

n = 3
G = 1

planets_mass = 1
dq0 = np.array([0.4662036850, 0.4662036850, -0.93240737, 0.4323657300, 0.4323657300, -0.86473146])

m = np.ones(n) * planets_mass
mm = m[:, None] * m[None, :]
q = np.array([-0.97000436, 0.97000436, 0, 0.24308753, -0.24308753, 0])
p = np.concatenate((dq0[:n] / m, dq0[n:] / m))

sizes = 1000 * np.arctan(m / 10)
soft = 0.01 * np.mean(np.abs(q)) # constant for softening

ANIMATION_LENGTH = 10
FRAMES_PER_SECOND = 60
FRAMES = ANIMATION_LENGTH * FRAMES_PER_SECOND
ANIMATION_INTERVAL = 1 / FRAMES_PER_SECOND
UPDATES_PER_FRAME = 1
DT = 1 / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)