import numpy as np

n = 100

M = 1000000
m = 1

max_T = 5
T = np.random.random(n-1) * max_T
r = ((T ** 2) * M / (4 * (np.pi ** 2))) ** (1 / 3)
v = (M / r) ** (1 / 2)
soft = 0.01 * np.mean(r) # epsilon for softening

m = np.array([M] + [m for _ in range(n-1)])
mm = m[:, None] * m[None, :]
q = np.concatenate((np.zeros((1, )), r, np.zeros((n, ))))
dq0 = np.concatenate((np.zeros((n + 1, )), v))
p = np.concatenate((dq0[:n] / m, dq0[n:] / m))

ANIMATION_LENGTH = 10
FRAMES_PER_SECOND = 60
FRAMES = ANIMATION_LENGTH * FRAMES_PER_SECOND
ANIMATION_INTERVAL = 1 / FRAMES_PER_SECOND
UPDATES_PER_FRAME = 1
DT = 1 / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)