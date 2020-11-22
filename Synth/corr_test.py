import matplotlib.pyplot as plt
import numpy as np
import datagen
import SWARMprocess


fs = 2
v = 7615

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

bobpos = np.array([150, 0])*v
bobpos = np.array([bobpos])
bobvel = np.array([694, 0])



bobvel = np.array([bobvel])

width = np.array([25])*v
A = np.array([1e5])
growth = np.array([0])

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
plt.show()

BA_shift = pro.timeshift_latitude(xposB, xposA, start = 0, stop = 500, shifts = 500)
BC_shift = pro.timeshift_latitude(xposB, xposC, start = 0, stop = 500, shifts = 500)
print("BA_shift = %g" % BA_shift)
print("BC_shift = %g" % BC_shift)

ind0 = 0
ind1 = 900
times = times[ind0:ind1]
dataB = dataB[ind0:ind1]
dataA = dataA[ind0 + BA_shift:ind1 + BA_shift]
dataC = dataC[ind0 + BC_shift:ind1 + BC_shift]

plt.plot(times, dataB)
plt.plot(times, dataA)
plt.plot(times, dataC)
plt.show()

diffB = dataB[1:] - dataB[:-1]
diffA = dataA[1:] - dataA[:-1]
diffC = dataC[1:] - dataC[:-1]

diffB = pro.meanie(diffB, 10)
diffA = pro.meanie(diffA, 10)
diffC = pro.meanie(diffC, 10)

plt.plot(diffB)
plt.plot(diffA)
plt.plot(diffC)
plt.show()

B_ind, BA_diff = pro.maxdiff(diffB, diffA)
B_ind2, BC_diff = pro.maxdiff(diffB, diffC)

s_ab = BA_shift
m_ab = s_ab + np.abs(BA_diff)

s_bc = BC_shift
m_bc = s_bc + np.abs(BC_diff)

print("BA velo = %G" % (v*(1 - s_ab/m_ab)))
print("BC velo = %g" % (v*(1 - s_bc/m_bc)))
