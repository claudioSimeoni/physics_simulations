import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

n = int(input("Insert the number of bodies: "))

M = 1000000
m = 1

max_T = 5
T = np.random.random(n-1) * max_T
r = ((T ** 2) * M / (4 * (np.pi ** 2))) ** (1 / 3)
v = (M / r) ** (1 / 2)
epsilon = 0.01 * np.mean(r) # epsilon for softening

m = np.array([M] + [m for _ in range(n-1)])
mm = m[:, None] * m[None, :]
q = np.concatenate((np.zeros((1, )), r, np.zeros((n, ))))
dq = np.concatenate((np.zeros((n + 1, )), v))
p = np.concatenate((dq[:n] / m, dq[n:] / m))


def grad_Hp():
    return np.concatenate((p[:n] / m, p[n:] / m))

def grad_Hq():

    x_diff = q[:n, None] - q[None, :n]
    y_diff = q[n:, None] - q[None, n:]

    module = (x_diff ** 2 + y_diff ** 2 + epsilon ** 2) ** 1.5
    np.fill_diagonal(module, np.inf)

    hq = np.concatenate((np.sum((mm * x_diff / module), axis=1), np.sum((mm * y_diff / module), axis=1)))

    return hq

dp = -grad_Hq()


# numerical_method
def leapfrog(): 
    global p, q, dp, dq
    pm = p + dt / 2 * dp
    q[:n] = q[:n] + dt * (p[:n] / m)
    q[n:] = q[n:] + dt * (p[n:] / m)
    dp = -grad_Hq()
    p = pm + dt / 2 * dp
    # dq = grad_Hp()


# animation constants 
animation_length = 10
frames_per_second = 60
frames = animation_length * frames_per_second
animation_interval = 1 / frames_per_second
updates_per_frame = 1
dt = 1 / (frames_per_second * updates_per_frame)


# plotting
xlim = np.max(r) * 100
ylim = xlim

fig, ax = plt.subplots()
fig.canvas.manager.full_screen_toggle()
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

ax.set_facecolor("black")
ax.set_xlim(-xlim, xlim)
ax.set_ylim(-ylim, ylim)

time = 0
sizes = 100 * np.arctan(m / 10)
bodies = ax.scatter(q[:n], q[n:], s=sizes, c="white")

def init(): 
    bodies.set_offsets(np.column_stack((q[:n], q[n:])))
    return bodies, 

def update(frame):
    global time
    for _ in range(updates_per_frame):
        leapfrog()
        time += dt

    bodies.set_offsets(np.column_stack((q[:n], q[n:])))
    return bodies, 

ani = animation.FuncAnimation(
    fig = fig, 
    func = update, 
    init_func = init, 
    frames = frames, 
    interval = animation_interval * 1000, # ms
    blit = True, 
    # repeat = False
)

plt.show()