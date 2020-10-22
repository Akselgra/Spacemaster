import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess as pro
from scipy.stats import pearsonr

pro = pro()

fs = 2
t = 100
f1 = 1/4
f2 = 1/3
f3 = 1/2
n = int(fs*t)
times = np.linspace(0, t, n)

samedat = np.sin(times*f1*2*np.pi)
data1 = samedat + np.sin(times*f2*2*np.pi)
data2 = samedat + np.sin(times*f3*2*np.pi)


fft1 = np.fft.fft(data1)/n
fft2 = np.fft.fft(data2)/n
fft1 = np.roll(fft1, int(n/2))
fft2 = np.roll(fft2, int(n/2))
freqs = np.linspace(-fs/2, fs/2, n)

plt.plot(freqs, np.abs(fft1))
plt.plot(freqs, np.abs(fft2))
plt.xlabel("Frequency [Hz]")
plt.ylabel("Power")
plt.legend(["1", "2"])
plt.title("Power spectrum of signal 1 and 2")
plt.show()

phase1 = np.arctan2(np.imag(fft1), np.real(fft2))
phase2 = np.arctan2(np.imag(fft2), np.real(fft2))

plt.plot(freqs, phase1)
plt.plot(freqs, phase2)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Phase")
plt.legend(["1", "2"])
plt.title("phase, arctan of imag/real")
plt.show()

cross_spec = pro.cross_spectrum(data1, data2, fs = fs)
cross_phase = np.arctan2(np.imag(cross_spec), np.real(cross_spec))

plt.plot(freqs, np.abs(cross_spec))
plt.xlabel("Frequency [Hz]")
plt.ylabel("Power")
plt.title("Power of cross spectrum")
plt.show()

plt.plot(freqs, cross_phase)
plt.xlabel("Frequency [Hz]")
plt.ylabel("Phase")
plt.title("phase, arctan of imag/real")
plt.show()

plt.plot(freqs, np.log10(np.abs(fft1)))
plt.plot(freqs, np.log10(np.abs(fft2)))
plt.plot(freqs, np.log10(np.abs(cross_spec)))
plt.xlabel("Frequency [Hz]")
plt.ylabel("log10(Power)")
plt.title("DTFT1, DTFT2, cross spectrum")
plt.legend(["1", "2", "CS"])
plt.show()
