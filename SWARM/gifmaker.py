import matplotlib.pyplot as plt
import numpy as np
from SWARMprocess import SWARMprocess
import imageio
import os
import glob
pro = SWARMprocess()

t = 100
fs = 2
n = int(t*fs)
times = np.linspace(0, t, n)
f0 = 0.1
f1 = 0.2
f2 = 0.3
frames = 100
fps = 15

data1 = (np.sin(times*f1*2*np.pi) + np.sin(times*f0*2*np.pi))/2
shifts = np.linspace(0.0, 2*np.pi, frames)

for shift in shifts:
    data2 = (np.sin(times*f2*2*np.pi + shift) + np.sin(times*f0*2*np.pi + shift))/2

    cross_spec = pro.cross_spectrum(data1, data2, fs =fs)
    cross_spec = np.roll(cross_spec, int(n/2))
    freqs = np.linspace(-fs/2, fs/2, n)
    phase = np.arctan2(np.imag(cross_spec), np.real(cross_spec))

    plt.plot(times, data1)
    plt.plot(times, data2)
    plt.xlabel("Time [s]")
    plt.ylabel("Amplitude")
    plt.legend(["stationary", "shifted"])
    plt.title("Signals, Shifted %2.2f" % shift)
    plt.axis([0, t, -1.1, 1.1])
    plt.savefig("signal_gif/signal_shift=%2.2f.png" % shift)
    plt.close()

    plt.plot(freqs, np.abs(cross_spec))
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Power")
    plt.title("CSD, shifted %2.2f" % shift)
    plt.axis([-fs/1.9, fs/1.9, np.min(np.abs(cross_spec)), np.max(np.abs(cross_spec))])
    plt.savefig("CSD_gif/CSD_shift=%2.2f.png" % shift)
    plt.close()

    plt.plot(freqs, phase)
    plt.axis([-fs/1.9, fs/1.9, -2*np.pi, 2*np.pi])
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Phase")
    plt.title("Phase, shifted %2.2f" % shift)
    plt.savefig("phase_gif/phase_shift=%2.2f.png" % shift)
    plt.close()



signal_dir = glob.glob("signal_gif/*.png")
CSD_dir = glob.glob("CSD_gif/*.png")
phase_dir = glob.glob("phase_gif/*.png")
signal_dir.sort()
CSD_dir.sort()
phase_dir.sort()


with imageio.get_writer("Figures/signals.gif", mode = "I", fps=fps) as writer:
    for signal in signal_dir:
        image = imageio.imread(signal)
        os.remove(signal)
        writer.append_data(image)

with imageio.get_writer("Figures/CSDs.gif", mode = "I", fps=fps) as writer:
    for CSD in CSD_dir:
        image = imageio.imread(CSD)
        os.remove(CSD)
        writer.append_data(image)

with imageio.get_writer("Figures/phases.gif", mode = "I", fps=fps) as writer:
    for phase in phase_dir:
        image = imageio.imread(phase)
        os.remove(phase)
        writer.append_data(image)
