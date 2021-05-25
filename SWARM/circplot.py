import numpy as np
import matplotlib.pyplot as plt

n = 100
r = np.linspace(0, 1, n)
r = np.zeros(n)+0.99
theta = np.linspace(0, 2*np.pi, n)

x0 = -0.25
y0 = -0.25
r0 = 0.5

circ_x = np.cos(theta)*r0 + x0
circ_y = np.sin(theta)*r0 + y0

circ_r = np.sqrt(circ_x**2 + circ_y**2)
circ_theta = np.arctan2(circ_y, circ_x)


fig, ax = plt.subplots(subplot_kw={"projection":"polar"})
ax.set_ylim(0, 1)

Thetas1 = 5*np.pi/4
Thetas2 = np.pi/4
Thetas = np.concatenate((np.zeros(int(n/2)) + Thetas1, np.zeros(int(n/2)) + Thetas2))
Rs = np.abs(np.linspace(-1, 1, n))
# ax.plot(theta, r)
ax.plot(circ_theta, circ_r)
ax.plot(Thetas, Rs)
plt.show()