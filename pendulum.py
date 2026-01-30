import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

# defining constants and numerical method's utility
xlim = 2
ylim = 1.5

l = ylim * 4 / 6
g = 9.81
xpivot = xlim / 2
ypivot = ylim * 5 / 6
initial_theta = np.pi / 8

animation_length = 10 # seconds
frames_per_second = 120
frames = animation_length * frames_per_second
animation_interval = 1 / frames_per_second
updates_per_frame = 5
dt = 1 / (frames_per_second * updates_per_frame)

def rk4(f, x, y, dx): 
    k1 = f(x, y)
    k2 = f(x + dx/2, y + dx/2 * k1)
    k3 = f(x + dx/2, y + dx/2 * k2)
    k4 = f(x + dx, y + dx * k3)
    return x + dx, y + dx/6 * (k1 + 2 * k2 + 2 * k3 + k4)

def f(x, y): 
    return np.array([y[1], - g / l * np.sin(y[0])])

# defining the plotting objects
fig, ax = plt.subplots()
ax.set_xlim(0, xlim)
ax.set_ylim(0, ylim)

# setting initial values of time, angle and coordinates of the ball 
t = 0
theta = np.array([initial_theta, 0]) # array of theta and theta'
xvalues = np.array([xpivot, xpivot + l * np.sin(initial_theta)])
yvalues = np.array([ypivot, ypivot - l * np.cos(initial_theta)])

pendulum, = ax.plot([], [], lw=2)
# ax.plot(np.linspace(0, xlim, xlim), np.ones((xlim, )) * yvalues[1])

def init(): 
    pendulum.set_data(xvalues, yvalues)
    return pendulum,

def update(frame): 
    global t, theta
    
    for _ in range(updates_per_frame):
        t, theta = rk4(f, t, theta, dt)

    xvalues[1] = xvalues[0] + l * np.sin(theta[0])
    yvalues[1] = yvalues[0] - l * np.cos(theta[0])

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