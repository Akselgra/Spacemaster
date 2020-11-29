import matplotlib.pyplot as plt
import numpy as np
import datagen
import SWARMprocess


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

bobvelx = np.random.random(n_bob)*0.1*v
bobvely = np.zeros_like(bobvelx)


bobpos = np.zeros((n_bob, 2))
bobvel = np.zeros((n_bob, 2))
for i in range(n_bob):
    bobpos[i][0] = bobposx[i]
    bobpos[i][1] = bobposy[i]
    bobvel[i][0] = bobvelx[i]
    bobvel[i][1] = bobvely[i]
    
width = np.zeros(n_bob) + 25*v
A = np.zeros(n_bob)+1e5
growth = np.zeros(n_bob)






n = int(fs*(t1 - t0))
times = np.linspace(t0, t1, n)
dataB = obj.satrun2(t0, t1, satpos1, satdir, bobpos, bobvel, width, A, growth)
dataA = obj.satrun2(t0, t1, satpos2, satdir, bobpos, bobvel, width, A, growth)
dataC = obj.satrun2(t0, t1, satpos3, satdir, bobpos, bobvel, width, A, growth)

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

plt.plot(times, dataB)
plt.plot(times, dataA)
plt.plot(times, dataC)
plt.xlabel("Time [s]")
plt.ylabel("Electron density")
plt.show()

