import numpy as np
import matplotlib.pyplot as plt
import SWARMprocess
pro = SWARMprocess.SWARMprocess()
N = 170000
n = 400
def bob(n):
    size = int(n/2)
    t = np.linspace(20, 160, 100)
    rands = []
    spatial = 0.2 + 0*(400-n)/400
    for i in range(len(t)):
        temporal = 0.002*t[i]
        rands.append(np.random.normal(scale = spatial + temporal, size = size))

    rands = np.array(rands)

    hists = []
    stds = []
    means = []
    bins = np.linspace(-1, 1, 50)

    for i in range(len(rands)):
        hist, bin = np.histogram(rands[i], bins = bins)
        hist = hist/np.sum(hist*(bin[1] - bin[0]))
        hists.append(hist)

    bins = (bins[1:] + bins[:-1])/2
    hists = np.array(hists)

    for i in range(len(rands)):
        temp_std, temp_mean = pro.std_mean(hists[i], bins)
        stds.append(temp_std)

    stds = np.array(stds)**2

    width = bins[1] - bins[0]

    p = np.polyfit(t, stds, deg = 1)
    a = p[0]
    b = p[1]

    return(stds, a, b)

ns = np.linspace(10, 400, 100)
a_s = []
b_s = []

for j in range(len(ns)):
    stds, a, b = bob(ns[j])
    a_s.append(a)
    b_s.append(b)

a_s = np.array(a_s)
b_s = np.array(b_s)

plt.figure(0)
plt.title("spatial")
plt.xlabel("n")
plt.ylabel("b")
plt.plot(ns, b_s)

plt.figure(1)
plt.title("temporal")
plt.xlabel("n")
plt.ylabel("a")
plt.plot(ns, a_s)

plt.figure(2)
plt.plot(ns, b_s/a_s)
plt.title("critical time")
plt.xlabel("n")
plt.show()
