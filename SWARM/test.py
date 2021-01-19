import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
pro = SWARMprocess()

std_0 = 0.25
x = np.random.normal(1, std_0, 10000)
mini = -2
maxy = 2

hist, bins = np.histogram(x, bins = np.linspace(mini, maxy, 50))

bins = (bins[1:] + bins[:-1])/2
width = bins[1] - bins[0]


hist = hist/np.sum(hist*width)


plt.bar(bins, hist, width = 0.9*width)

mean = np.sum(hist*bins)/np.sum(hist)
#mean = np.average(bins, weights = hist)

std = np.sqrt(np.average((bins-mean)**2, weights = hist))


thing = (bins-mean)**2
std = np.sqrt(np.sum(thing*hist)/np.sum(hist))



std, mean = pro.std_mean(hist, bins)
y = np.linspace(mini, maxy, 1000)
gauss = pro.gauss_curve(y, mean, std)
plt.plot(y, gauss, "r")
plt.show()

plt.plot(bins, bins*hist, "o")
plt.show()