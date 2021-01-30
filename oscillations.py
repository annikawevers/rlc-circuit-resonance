import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

'''
The following program generates resonance plots of an RLC resonant circuit
for underdamped, critically damped, and overdamped cases. These plots show
current vs frequency. In addition, the time evolution of the current for each
of the three cases is also plotted.
'''


L = input('Enter value for L: ')
C = input('Enter value for C: ')
L = float(L)
C = float(C)

# Define initial conditions
w0 = 1/(np.sqrt(L*C))
ta = np.linspace(0, 100, 3000 + 1)
Run = 0.2*2*np.sqrt(L/C)
Rcr = 2*np.sqrt(L/C)
Rov = 10*np.sqrt(L/C)
y0a = [0, 0]
Rarray = np.array([10*np.sqrt(L/C), 2*np.sqrt(L/C), 0.2*2*np.sqrt(L/C)])
wvals = np.linspace(0.1*w0, 2*w0, 50 + 1)
wplot = np.linspace(0.1*w0/(2*np.pi), 2*w0/(2*np.pi), 50+1)


# Differential equation
def f(ya, ta, L, C, R, w):
    x = ya[0]
    v = ya[1]
    dvdt = (1/L)*(np.cos(w*ta)) - (1/(L*C))*x - (R/L)*v
    return np.array([v, dvdt])


# Underdamped case:
max = []
for k in wvals:
    sol = odeint(f, y0a, ta, args=(L, C, Run, k))
    Ivals = sol[200:, 1]
    Imax = np.amax(Ivals)
    max.append(Imax)
y1 = 2*np.array(max)
plt.subplot(311)
plt.plot(wplot, y1)
plt.title('Max current vs frequency (Underdamped)')
plt.xlabel('frequency (Hz)')
plt.ylabel('current (A)')

# Critically damped case:
max2 = []
for k in wvals:
    sol2 = odeint(f, y0a, ta, args=(L, C, Rcr, k))
    Ivals2 = sol2[200:, 1]
    Imax2 = np.amax(Ivals2)
    max2.append(Imax2)
y2 = 2*np.array(max2)
plt.subplot(312)
plt.plot(wplot, y2)
plt.title('Max current vs frequency (Critically Damped)')
plt.xlabel('frequency (Hz)')
plt.ylabel('current (A)')

# Overdamped case:
max3 = []
for k in wvals:
    sol3 = odeint(f, y0a, ta, args=(L, C, Rov, k))
    Ivals3 = sol3[200:, 1]
    Imax3 = np.amax(Ivals3)
    max3.append(Imax3)
y3 = 2*np.array(max3)
plt.subplot(313)
plt.plot(wplot, y3)
plt.title('Max current vs frequency (Overdamped)')
plt.xlabel('frequency (Hz)')
plt.ylabel('current (A)')

plt.tight_layout()
plt.savefig('resonance.pdf')
plt.clf()


# Time evolution of I(t)
ta = np.linspace(0, 75, 375+1)  # new time array for new plots

plt.subplot(311)  # underdamped case
sol = odeint(f, y0a, ta, args=(L, C, Run, w0))
y1 = sol[:, 1]
plt.plot(ta, y1)
plt.title('Evolution of current (underdamped)')
plt.xlabel('time (s)')
plt.ylabel('current (A)')

plt.subplot(312)  # critically damped case
sol2 = odeint(f, y0a, ta, args=(L, C, Rcr, w0))
y2 = sol2[:, 1]
plt.plot(ta, y2)
plt.title('Evolution of current (critically damped)')
plt.xlabel('time (s)')
plt.ylabel('current (A)')

plt.subplot(313)  # overdamped case
sol3 = odeint(f, y0a, ta, args=(L, C, Rov, w0))
y3 = sol3[:, 1]
plt.plot(ta, y3)
plt.title('Evolution of current (overdamped)')
plt.xlabel('time (s)')
plt.ylabel('current (A)')

plt.tight_layout()  # spaces plots evenly
plt.savefig('transients.pdf')
