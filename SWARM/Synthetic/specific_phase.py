import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess as pro
from scipy.stats import pearsonr

pro = pro()

fs = 2
t = 100
n = int(fs*t)
f0 = 0.2
f1 = 0.3
f2 = 0.4
times = np.linspace(0, t, n)

n_shifts = 100
maxshift = 9.99
shifts = np.linspace(0, maxshift, n_shifts)


data1 = (np.sin(times*f0*2*np.pi) + np.sin(times*f1*2*np.pi))/2

phases = []

for shift in shifts:
    times2 = np.linspace(shift, t + shift, n)
    data2 = (np.sin(times2*f0*2*np.pi) + np.sin(times2*f2*2*np.pi))/2
    cross_spec = pro.cross_spectrum(data1, data2)[:int(n/2)]
    phase = np.arctan2(np.imag(cross_spec), np.real(cross_spec)) + np.pi
    freqs = np.linspace(0, fs/2, int(n/2))
    ind = int(np.round(f0/fs*2*((n/2)-1)))
    phases.append(phase[ind])


plt.plot(shifts, phases)
plt.xlabel("Time shifted")
plt.ylabel("phase of cross spectrum")
plt.title("Cross spectrum phase at frequency = %g" % f0)
plt.show()
