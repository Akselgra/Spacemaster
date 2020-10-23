import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
pro = SWARMprocess()

fs = 2
t = 100
n = int(fs*t)
f = 0.1
f1 = 0.2
f2 = 0.3
times = np.linspace(0, t, n)
shift = 3

base = np.sin(times*f*np.pi*2)

data1 = np.sin(times*f1*np.pi*2) + base
data2 = np.sin(times*f2*np.pi*2) + base

plt.plot(times, data1)
plt.plot(times, data2)
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.title("Sine functions")
plt.legend(["f = %g Hz" % f1, "f = %g Hz" % f2])
plt.show()

corrvec = pro.correlator2(data1, data2)
shiftvec = np.linspace(-t/2, t/2, n)

plt.plot(shiftvec, corrvec)
plt.plot([0,0],[-1.1,1.1], "k")
plt.xlabel("Time shifted")
plt.ylabel("Pearson R correlation coefficient")
plt.title("cross correlations, before shift")
plt.show()

cross_spec = pro.cross_spectrum(data1, data2, fs = fs)
cross_spec = np.roll(cross_spec, int(n/2))
phase = np.arctan2(np.imag(cross_spec), np.real(cross_spec))
freqs = np.linspace(-fs/2, fs/2, n)

plt.plot(freqs, np.abs(cross_spec))
plt.xlabel("Frequency [Hz]")
plt.ylabel("Power density")
plt.title("Cross spectrum density")
plt.show()

plt.plot(freqs, phase)
plt.xlabel("Frequency [Hz]")
plt.ylabel("phase")
plt.title("Cross spectrum phase")
plt.show()

data2 = np.roll(data2, int(shift*fs))

corrvec = pro.correlator2(data1, data2)

plt.plot(times, data1)
plt.plot(times, data2)
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")
plt.title("sines, shifted")
plt.legend(["f = %g Hz" % f1, "f = %g Hz" % f2])
plt.show()


plt.plot(shiftvec, corrvec)
plt.plot([0,0],[-1.1, 1.1], "k")
plt.xlabel("Time shifted [s]")
plt.ylabel("Pearson R correlation coefficient")
plt.title("Cross correlations, after shift")
plt.show()

cross_spec = pro.cross_spectrum(data1, data2, fs = fs)
cross_spec = np.roll(cross_spec, int(n/2))
phase = np.arctan2(np.imag(cross_spec), np.real(cross_spec))
freqs = np.linspace(-fs/2, fs/2, n)

plt.plot(freqs, np.abs(cross_spec))
plt.xlabel("Frequency [Hz]")
plt.ylabel("Power density")
plt.title("Cross spectrum density")
plt.show()

plt.plot(freqs, phase)
plt.xlabel("Frequency [Hz]")
plt.ylabel("phase")
plt.title("Cross spectrum phase")
plt.show()
