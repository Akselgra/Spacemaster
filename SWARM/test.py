import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess as pro
from scipy.stats import pearsonr

pro = pro()

t = 3
fs = 1000
n = int(fs*t)

zeros = np.zeros(int(n/3))
window = np.concatenate((zeros, zeros+1, zeros))

plt.plot(window)
plt.show()

corr = pro.correlator2(window, window)


corr2 = np.correlate(window, window, "full")

plt.plot(corr)
plt.plot(corr2/corr2.max())
plt.show()
