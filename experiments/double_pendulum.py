import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation

# defining constants and numerical method's utility
xlim = 7
ylim = 4

l1 = ylim * 2 / 6
l2 = ylim * 2 / 6
g = 9.81
m1 = 20
m2 = 100
xpivot = xlim / 2
ypivot = ylim * 5 / 6
initial_theta1 = np.pi / 8
initial_theta2 = np.pi

animation_length = 10 # seconds
frames_per_second = 60
frames = animation_length * frames_per_second
animation_interval = 1 / frames_per_second
updates_per_frame = 1
dt = 1 / (frames_per_second * updates_per_frame)

def rk4(f, x, y, dx): 
    k1 = f(x, y)
    k2 = f(x + dx/2, y + dx/2 * k1)
    k3 = f(x + dx/2, y + dx/2 * k2)
    k4 = f(x + dx, y + dx * k3)
    return x + dx, y + dx/6 * (k1 + 2 * k2 + 2 * k3 + k4)

# th1, th2, th1', th2'
# def f(x, y): 
#     return np.array([y[2], 
#                      y[3], 
#                      g * (m2 * y[1] - (m1 + m2) * y[0]) / (m1 * l1), 
#                      g * (m1 + m2) * (y[0] - y[1]) / (m1 * l2)])
def f(x, y):
    th1, th2, w1, w2 = y
    delta = th2 - th1
    den1 = (m1 + m2) * l1 - m2 * l1 * np.cos(delta) * np.cos(delta)
    den2 = (l2/l1) * den1

    dth1 = w1
    dth2 = w2
    dw1 = (m2*l1*w1*w1*np.sin(delta)*np.cos(delta) + m2*g*np.sin(th2)*np.cos(delta) + m2*l2*w2*w2*np.sin(delta) - (m1+m2)*g*np.sin(th1)) / den1
    dw2 = (-m2*l2*w2*w2*np.sin(delta)*np.cos(delta) + (m1+m2)*(g*np.sin(th1)*np.cos(delta) - l1*w1*w1*np.sin(delta) - g*np.sin(th2))) / den2

    return np.array([dth1, dth2, dw1, dw2])

# defining the plotting objects
fig, ax = plt.subplots()
ax.set_xlim(0, xlim)
ax.set_ylim(0, ylim)

# setting initial values of time, angle and coordinates of the ball 
t = 0
theta = np.array([initial_theta1, initial_theta2, 0, 0]) # array of theta and theta'
x1 = xpivot + l1 * np.sin(initial_theta1)
y1 = ypivot - l1 * np.cos(initial_theta1)
xvalues = np.array([xpivot, x1, x1 + l2 * np.sin(initial_theta2)])
yvalues = np.array([ypivot, y1, y1 - l2 * np.cos(initial_theta2)])

pendulum, = ax.plot([], [], lw=2, color="black", solid_capstyle="round")
mass1, = ax.plot([], [], 'o', color="blue", markersize=40 * np.sqrt(1))
mass2, = ax.plot([], [], 'o', color="green", markersize=40 * np.sqrt(m2 / m1))


def init(): 
    pendulum.set_data(xvalues, yvalues)
    mass1.set_data(xvalues[1], yvalues[1])
    mass2.set_data(xvalues[2], yvalues[2])
    return pendulum, mass1, mass2

def update(frame): 
    global t, theta
    
    for _ in range(updates_per_frame):
        t, theta = rk4(f, t, theta, dt)

    xvalues[1] = xvalues[0] + l1 * np.sin(theta[0])
    yvalues[1] = yvalues[0] - l1 * np.cos(theta[0])
    xvalues[2] = xvalues[1] + l2 * np.sin(theta[1])
    yvalues[2] = yvalues[1] - l2 * np.cos(theta[1])

    pendulum.set_data(xvalues, yvalues)
    mass1.set_data(xvalues[1], yvalues[1])
    mass2.set_data(xvalues[2], yvalues[2])
    return pendulum, mass1, mass2

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