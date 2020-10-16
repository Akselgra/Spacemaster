import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess as pro
from scipy.stats import norm

pro = pro()

n = 50000
f1 = 2*np.pi
f2 = 2*np.pi
t = 10
fs = n/t
phase1 = np.pi
phase2 = 0
A = 1

data1 = A*np.sin(np.linspace(phase1, f1*t + phase1, n)) + np.random.random(n)
data2 = A*np.sin(np.linspace(phase2, f2*t + phase2, n)) + np.random.random(n)

data1 = pro.meanie(data1, 200)
data2 = pro.meanie(data2, 200)
plt.plot(data1)
plt.plot(data2)
plt.show()

corr = np.correlate(data1, data2, "full")
indices = np.arange(len(corr))
times = np.linspace(-t, t, len(corr))

plt.plot(times, corr)
plt.show()
