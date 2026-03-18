import numpy as np

n = 100
G = 1

star_mass = 1000000
planets_mass = 1

max_T = 5
T = np.random.random(n - 1) * max_T
r = ((T ** 2) * star_mass * G / (4 * (np.pi ** 2))) ** (1 / 3)
v = (star_mass * G / r) ** (1 / 2)
dq0 = np.concatenate((np.zeros((n + 1, )), v))


m = np.array([star_mass] + [planets_mass for _ in range(n-1)])
mm = m[:, None] * m[None, :]
q = np.concatenate((np.zeros((1, )), r, np.zeros((n, ))))
p = np.concatenate((dq0[:n] / m, dq0[n:] / m))

sizes = 100 * np.arctan(m / 10)
soft = 0.01 * np.mean(np.abs(q))

ANIMATION_LENGTH = 10
FRAMES_PER_SECOND = 60
FRAMES = ANIMATION_LENGTH * FRAMES_PER_SECOND
ANIMATION_INTERVAL = 1 / FRAMES_PER_SECOND
UPDATES_PER_FRAME = 1
DT = 1 / (FRAMES_PER_SECOND * UPDATES_PER_FRAME)