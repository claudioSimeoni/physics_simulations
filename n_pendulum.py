from sympy import * 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# symbolic expressions definition and reading initial values
n = int(input())

t = symbols("t")
th = symbols(f"th0:{n}", cls=Function)
l = symbols(f"l0:{n}")
m = symbols(f"m0:{n}")
g = symbols("g")

xlim = 0.25
ylim = 0.30
xpivot = 0
ypivot = 0

print("Insert l values: ")
lvals = [xlim / 6 for _ in range(n)]
print("Insert m values: ")
mvals = [1 for _ in range(n)]
print("Insert th values: ")
thvals = [1 + 0.1 * i for i in range(n)]
gval = 9.81

x = [sum(l[k] * sin(th[k](t)) for k in range(i + 1)) for i in range(n)]
y = [-sum(l[k] * cos(th[k](t)) for k in range(i + 1)) for i in range(n)]

vx = [diff(x[i], t) for i in range(n)]
vy = [diff(y[i], t) for i in range(n)]

T = sum(m[i] * (vx[i] ** 2 + vy[i] ** 2) / 2 for i in range(n))
U = g * sum(m[i] * y[i] for i in range(n))
L = T - U 

eq = [diff(diff(L, diff(th[i](t), t)), t) - diff(L, th[i](t)) for i in range(n)]

# replacing in the equation
th_sym   = symbols(f'th_sym0:{n}')
thd_sym  = symbols(f'thd_sym0:{n}')
thdd_sym = symbols(f'thdd_sym0:{n}')

rep = {g: gval}
for j in range(n):
    rep[th[j](t)] = th_sym[j]
    rep[Derivative(th[j](t), t)] = thd_sym[j]
    rep[Derivative(th[j](t), t, 2)] = thdd_sym[j]
    rep[m[j]] = mvals[j]
    rep[l[j]] = lvals[j]

eq_subs = [eq[i].xreplace(rep) for i in range(n)]
# eq_resp_dd = list(linsolve(eq_subs, [thdd_sym[j] for j in range(n)]))[0]
# eq_function = lambdify([thd_sym, th_sym], eq_resp_dd, "numpy")

# constants 
animation_length = 10
frames_per_second = 60
frames = animation_length * frames_per_second
animation_interval = 1 / frames_per_second
updates_per_frame = 1
dt = 1 / (frames_per_second * updates_per_frame)

# numerical_method
def rk4(f, x, y, dx): 
    k1 = f(x, y)
    k2 = f(x + dx/2, y + dx/2 * k1)
    k3 = f(x + dx/2, y + dx/2 * k2)
    k4 = f(x + dx, y + dx * k3)
    return x + dx, y + dx/6 * (k1 + 2 * k2 + 2 * k3 + k4)

def f(time, pos):
    subs_list = (
        {thd_sym[j]: pos[j+n] for j in range(n)} |
        {th_sym[j]: pos[j] for j in range(n)}
    )
  
    eq_subs_curr = [eq_subs[i].xreplace(subs_list) for i in range(n)]
    new_th = pos[n:]
    new_thd = list(linsolve(eq_subs_curr, [thdd_sym[j] for j in range(n)]))[0]
    return np.concatenate((np.array(new_th), np.array(new_thd)))

# plotting
fig, ax = plt.subplots()
ax.set_xlim(-xlim, xlim)
ax.set_ylim(-ylim, 0)

# setting initial values of time and angle
time = 0
pos = np.array(thvals + [0 for _ in range(n)])
xvalues = np.array([xpivot] + [sum(lvals[k] * np.sin(thvals[k]) for k in range(i + 1)) for i in range(n)])
yvalues = np.array([ypivot] + [-sum(lvals[k] * np.cos(thvals[k]) for k in range(i + 1)) for i in range(n)])

pendulum, = ax.plot([], [], lw=2, color="black", solid_capstyle="round")

def init(): 
    pendulum.set_data(xvalues, yvalues)
    return pendulum, 

def update(frame): 
    global time, pos
    
    for _ in range(updates_per_frame):
        time, pos = rk4(f, time, pos, dt)
        pos = np.array([float(x) for x in pos])

    xvalues = np.array([xpivot] + [sum(lvals[k] * np.sin(pos[k]) for k in range(i + 1)) for i in range(n)])
    yvalues = np.array([ypivot] + [-sum(lvals[k] * np.cos(pos[k]) for k in range(i + 1)) for i in range(n)])

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