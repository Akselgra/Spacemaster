import matplotlib.pyplot as plt
import numpy as np
import datagen
import SWARMprocess
import scipy.signal.windows as windows


pro = SWARMprocess.SWARMprocess()
fs = 2
v = 7615
dx = v/fs #distance between data points

obj = datagen.SynthDat(fs = fs, v = v)
pro = SWARMprocess.SWARMprocess()
#t0, t1, satpos, satdir, bobpos, bobvel, width, A, growth

angle = 0 # angle of satellite direction in degrees
angle = np.deg2rad(angle)
satdir = np.array([np.cos(angle), np.sin(angle)])

t0 = 0
t1 = 1000

satpos = np.array([0, 0])*v




bobpos = np.array([100, 0])*v
bobpos = np.array([bobpos])
bobvel = np.array([0, 0])
bobvel = np.array([bobvel])
    
width = np.array([10])*v
A = np.array([1e5])
growth = np.array([0])

n = int((t1-t0)*fs)
data1 = obj.satrun2(t0, t1, satpos, satdir, bobpos, bobvel, width, A, growth, noise = False)
times = np.linspace(t0, t1, n)



freqs = np.linspace(0, fs/2, int(n/2))


width2 = np.array([20])*v

data2 = obj.satrun2(t0, t1, satpos, satdir, bobpos, bobvel, width2, A, growth, noise = False)

width3 = np.array([30])*v
data3 = obj.satrun2(t0, t1, satpos, satdir, bobpos, bobvel, width3, A, growth, noise = False)

plt.plot(times, data1)
plt.plot(times, data2)
plt.plot(times, data3)
plt.xlabel("Time [s]")
plt.ylabel("Electron density")
plt.legend(["Width = 10 s*v", "width = 20s*v", "width = 30 s*v"])
plt.show()

fft1 = np.fft.fft(data1)[:int(n/2)]
fft2 = np.fft.fft(data2)[:int(n/2)]
fft3 = np.fft.fft(data3)[:int(n/2)]

plt.figure(1)
plt.plot(freqs, np.abs(fft1))
plt.plot(freqs, np.abs(fft2))
plt.plot(freqs, np.abs(fft3))
plt.xlabel("Frequency [Hz]")
plt.ylabel("Power spectrum density")
plt.legend(["Width = 10 s*v", "width = 20s*v", "width = 30 s*v"])
plt.title("PSD")
#plt.show()

CSD_12 = pro.cross_spectrum(data1, data2)
CSD_13 = pro.cross_spectrum(data1, data3)
CSD_23 = pro.cross_spectrum(data2, data3)
plt.figure(2)
plt.plot(freqs, np.abs(CSD_12[:int(n/2)]))
plt.plot(freqs, np.abs(CSD_13[:int(n/2)]))
plt.plot(freqs, np.abs(CSD_23[:int(n/2)]))
plt.xlabel("Frequency [Hz]")
plt.ylabel("Cross spectrum density")
plt.legend(["10-20", "10-30", "20-30"])
plt.title("CSD")
plt.show()
