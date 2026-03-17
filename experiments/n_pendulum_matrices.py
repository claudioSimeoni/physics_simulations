from sympy import * 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# symbolic expressions definition and reading initial values
n = int(input())

th = symbols(f'th0:{n}')
thd = symbols(f'thd0:{n}')

xlim = 4
ylim = 3
xpivot = 0
ypivot = 0

print("Insert l values: ")
l = np.array([ylim / n * 0.8 for _ in range(n)], dtype=np.float64)
print("Insert m values: ")
m = np.array([1 for _ in range(n)], dtype=np.float64)
print("Insert th values: ")
th0 = np.array([1 + 0.1 * i for i in range(n)], dtype=np.float64)
g = 9.81


C = Matrix([sum(
                sum(m[k] * l[i] * l[j] * sin(th[i] - th[j]) * (thd[j] ** 2) 
                    for k in range(max(i, j), n))
                for j in range(n))
            for i in range(n)])

G = Matrix([g * l[i] * sin(th[i]) * 
                sum(m[k] for k in range(i, n)) 
            for i in range(n)])
M = Matrix([[l[i] * l[j] * 
                sum(m[k] * cos(th[i] - th[j]) for k in range(max(i, j), n)) 
                for j in range(n)] 
            for i in range(n)])

M_func = lambdify((th, thd), M, "numpy")
b_func = lambdify((th, thd), -C - G, "numpy")


# constants 
animation_length = 10
frames_per_second = 120
frames = animation_length * frames_per_second
animation_interval = 1 / frames_per_second
updates_per_frame = 10
dt = 1 / (frames_per_second * updates_per_frame)

# numerical_method
def rk4(f, x, y, dx): 
    k1 = f(x, y)
    k2 = f(x + dx/2, y + dx/2 * k1)
    k3 = f(x + dx/2, y + dx/2 * k2)
    k4 = f(x + dx, y + dx * k3)
    return x + dx, y + dx/6 * (k1 + 2 * k2 + 2 * k3 + k4)

def f(t, pos):
    M_t = M_func(pos[:n], pos[n:])
    b_t = b_func(pos[:n], pos[n:])
    new_pos = pos[n:]
    new_posd = np.linalg.solve(M_t, b_t).reshape(n,)
    return np.concatenate([new_pos, new_posd])

# plotting
fig, ax = plt.subplots()
ax.set_xlim(-xlim, xlim)
ax.set_ylim(-ylim, ylim / 2)

# setting initial values of time and angle
time = 0
pos = np.array(list(th0) + [0 for _ in range(n)])
xvalues = np.array([xpivot] + [sum(l[k] * np.sin(th0[k]) for k in range(i + 1)) for i in range(n)])
yvalues = np.array([ypivot] + [-sum(l[k] * np.cos(th0[k]) for k in range(i + 1)) for i in range(n)])

pendulum, = ax.plot([], [], lw=2, color="black", solid_capstyle="round")

def init(): 
    pendulum.set_data(xvalues, yvalues)
    return pendulum, 

def update(frame): 
    global time, pos
    
    for _ in range(updates_per_frame):
        time, pos = rk4(f, time, pos, dt)

    xvalues = np.array([xpivot] + [sum(l[k] * np.sin(pos[k]) for k in range(i + 1)) for i in range(n)])
    yvalues = np.array([ypivot] + [-sum(l[k] * np.cos(pos[k]) for k in range(i + 1)) for i in range(n)])

    pendulum.set_data(xvalues, yvalues)
    return pendulum, 

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