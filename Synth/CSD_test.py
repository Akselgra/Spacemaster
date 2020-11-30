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

satpos1 = np.array([0, 0])*v
satpos2 = np.array([-45, 0])*v
satpos3 = np.array([-123, 0])*v

n_bob = 5
bobposx = np.linspace(50, 450, n_bob)*v
bobposy = np.zeros_like(bobposx)

bobvelx = np.random.random(n_bob)*0*v
bobvely = np.zeros_like(bobvelx)


bobpos = np.zeros((n_bob, 2))
bobvel = np.zeros((n_bob, 2))
for i in range(n_bob):
    bobpos[i][0] = bobposx[i]
    bobpos[i][1] = bobposy[i]
    bobvel[i][0] = bobvelx[i]
    bobvel[i][1] = bobvely[i]
    
width = np.zeros(n_bob) + 10*v
A = np.zeros(n_bob)+1e5*0.8 + np.random.random(n_bob)*1e5*0.2
growth = np.zeros(n_bob) + 75






n = int(fs*(t1 - t0))
times = np.linspace(t0, t1, n)
dataB = obj.satrun2(t0, t1, satpos1, satdir, bobpos, bobvel, width, A, growth, noise = True)
dataA = obj.satrun2(t0, t1, satpos2, satdir, bobpos, bobvel, width, A, growth, noise = True)
dataC = obj.satrun2(t0, t1, satpos3, satdir, bobpos, bobvel, width, A, growth, noise = True)

xposB = times*v*satdir[0] + satpos1[0]
xposA = times*v*satdir[0] + satpos2[0]
xposC = times*v*satdir[0] + satpos3[0]


plt.plot(times, dataB)
plt.plot(times, dataA)
plt.plot(times, dataC)
plt.xlabel("Time [s]")
plt.ylabel("Electron density")
plt.show()

BA_shift = pro.timeshift_latitude(xposB, xposA, start = 0, stop = 500, shifts = 500)
BC_shift = pro.timeshift_latitude(xposB, xposC, start = 0, stop = 500, shifts = 500)
print("BA_shift = %g" % BA_shift)
print("BC_shift = %g" % BC_shift)

ind0 = 0
ind1 = 1200
times = times[ind0:ind1]
dataB = dataB[ind0:ind1]
dataA = dataA[ind0 + BA_shift:ind1 + BA_shift]
dataC = dataC[ind0 + BC_shift:ind1 + BC_shift]

plt.figure(1)
plt.plot(times, dataB)
plt.plot(times, dataA)
plt.plot(times, dataC)
plt.xlabel("Time [s]")
plt.ylabel("Electron density")
plt.axis([0, 600, 0, 1.25e5])
#plt.show()

#Creating CSDs

t0 = 50
t1 = 600
n = 25

times = np.linspace(t0, t1, n)
indices = np.zeros_like(times)
for i in range(len(indices)):
    indices[i] = int(times[i]*fs)

    
n_freq = int(indices[1] - indices[0])
freqs = np.linspace(0, fs/2, n_freq)
BA_CSDs = []
BC_CSDs = []
AC_CSDs = []

for i in range(len(indices)-1):
    ind0 = int(indices[i] - n_freq)
    ind1 = int(indices[i+1])
    
    n_window = ind1 - ind0
    window = windows.general_gaussian(n_window, 1, sig = n_window/8)
    curr_dataA = dataA[ind0:ind1]
    curr_dataB = dataB[ind0:ind1]
    curr_dataC = dataC[ind0:ind1]
    
    curr_dataA = curr_dataA - np.mean(curr_dataA)
    curr_dataB = curr_dataB - np.mean(curr_dataB)
    curr_dataC = curr_dataC - np.mean(curr_dataC)
    
    curr_dataA = curr_dataA*window
    curr_dataB = curr_dataB*window
    curr_dataC = curr_dataC*window
    
    BA_cross_spec = pro.cross_spectrum(curr_dataB, curr_dataA)[:n_freq]
    BC_cross_spec = pro.cross_spectrum(curr_dataB, curr_dataC)[:n_freq]
    AC_cross_spec = pro.cross_spectrum(curr_dataA, curr_dataC)[:n_freq]

    BA_CSDs.append(np.abs(BA_cross_spec))
    BC_CSDs.append(np.abs(BC_cross_spec))
    AC_CSDs.append(np.abs(AC_cross_spec))
    

BA_CSDs = np.array(BA_CSDs) + 1e-16
BC_CSDs = np.array(BC_CSDs) + 1e-16
AC_CSDs = np.array(AC_CSDs) + 1e-16


Freqs, Times = np.meshgrid(freqs, times[:-1])


loggy = np.log10(BA_CSDs)
#loggy = BA_CSDs
plt.figure(2)
plt.pcolormesh(Times, Freqs, loggy, cmap = "magma")
plt.colorbar()
plt.xlabel("Time [s]")
plt.ylabel("Frequency [Hz]")
plt.show()
    
