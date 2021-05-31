import numpy as np
import matplotlib.pyplot as plt

def square(x, mid, width, val = 100):
    """
    defines a square pulse signal
    """
    temp_arr = np.copy(x)
    inds1 = np.nonzero(np.abs(x-mid) <= width)[0]
    inds2 = np.nonzero(np.abs(x-mid) > width)[0]
    temp_arr[inds1] = val
    temp_arr[inds2] = 0
    return temp_arr


mid = 2
width = 1
t = 10
fs = 100
n = int(fs*t)
df = 1/t
print(df)
val = 10000

xs = np.linspace(0, t, n)

sig = square(xs, mid, width, val)
sig = sig + np.random.normal(size = len(sig), scale = 0.01*val)

plt.plot(xs, sig)
plt.show()

fft = np.fft.fft(sig)[:int(n/2)]/n*2
PSD = np.abs(fft)**2
freqs = np.arange(int(n/2))*df

freqs = freqs[1:]
PSD = PSD[1:]

PSD += 1e-16


freqs = np.log10(freqs)
PSD = np.log10(PSD)
plt.plot(freqs, PSD, ".")


a, b = np.polyfit(freqs, PSD, deg = 1)
plt.plot(freqs, a*freqs + b)
ticks = np.array([freqs[0], freqs[5], freqs[20], freqs[300], freqs[-1]])
strs = []
for i in range(len(ticks)):
    strs.append("%g" % (10**ticks[i]))
plt.xticks(ticks = ticks, labels = strs)
plt.show()


print(a)
print(5/3)


