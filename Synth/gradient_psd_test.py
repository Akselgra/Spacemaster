import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess as pro
pro = pro()


def gradient(x, start, stop, val = 100):
    """
    input:
        x - current position [float]
        start - position where the gradient starts [float]
        stop - position where gradient stops [float]
        val - value at stop [float]
    output:
        value if start < x < stop:
        0 else
    """

    X = (x - start)
    if start < x < stop:
        return(X/(stop - start)*val)
    else:
        return(0)

def triangle(x, mid, width, val = 100):
    """
    input:
        x - current position [float]
        mid - position of triangle [float]
        width - width of triangle [float]
        val - max value of triangle [float]
    returns
        triangle value at x
    """

    X = np.abs(np.abs(x - mid) - width)
    if np.abs(x - mid) < width:
        return(X/width*val)
    else:
        return 0



t = 1000
fs = 2
n = int(fs*t)
xs = np.linspace(0, t, n)
mid = 300
width = 50

#xs = np.linspace(0, 1000, 5000)
ys = np.zeros_like(xs)

mids = [100, 300, 500, 700]
widths = [3, 5, 10, 20]
for j in range(len(mids)):
    for i in range(len(xs)):
        ys[i] += triangle(xs[i], mids[j], widths[j])

# omegas = np.linspace(0, t, n)*0.5*2*np.pi
# sin = np.sin(omegas)*10
# ys += sin
ys += np.random.normal(size = len(ys))

Freqs, Times, ffts = pro.fft_time(ys, 100, fs)

ffts = np.abs(ffts) + 1e-16
#ffts = np.log10(ffts)

# plt.figure(1)
# plt.plot(xs, ys)
# plt.grid(b = True)
#
# plt.figure(2)
# plt.pcolormesh(Times, Freqs, np.log10(ffts), cmap = "magma")
# plt.colorbar()


freqs = Freqs[0, :]
df = freqs[1] - freqs[0]
times = Times[:, 0]
# fourier_ints = np.zeros_like(times)
# #print(len(fourier_int))
# print(len(Times[:, 0]))
# print(len(ffts[:, 0]))
# for i in range(len(times)):
#     fourier_ints[i] = np.sum(ffts[:int(len(freqs)/8), i])*df


temp_ffts1 = ffts[:int(len(times)), :]
temp_ffts2 = ffts[:int(len(times)), :5]
fourier_int1 = np.sum(temp_ffts1, axis = 1)*df
fourier_int2 = np.sum(temp_ffts2, axis = 1)*df

print(np.shape(temp_ffts1))
print(np.shape(Times[:, 0]))


# plt.figure(3)
# plt.plot(Times[:, 0], fourier_int1)
# plt.grid(b = True)
#
# plt.figure(4)
# plt.plot(Times[:, 0], fourier_int2)
# plt.grid(b = True)
# plt.show()

plt.figure(1)
plt.plot(xs, ys/np.max(ys))
plt.plot(Times[:, 0], fourier_int2/np.max(fourier_int1))

Freqs, Times, ffts = pro.fft_time(ys, width, fs)
plt.figure(2)
plt.pcolormesh(Times, Freqs, np.log10(np.abs(ffts)), cmap = "magma")
plt.colorbar()

minfreq = 0.25
maxfreq = 1

times, fft_ints = pro.fft_time_integral(ys, width, fs, minfreq, maxfreq)
plt.figure(3)
plt.plot(times, fft_ints)
plt.show()
