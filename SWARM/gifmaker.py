import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation




t = 100
fs = 2
n = int(t*fs)
times = np.linspace(0, t, n)
frames = 200

fig = plt.figure()
ax = plt.axes(xlim=(0, t), ylim = (-1.1, 1.1))
thing, = ax.plot([], [])


f1 = 0.1
f2 = 0.1


shifts = np.linspace(0, 4*np.pi, 4)

def init():
    data1 = np.sin(times*2*np.pi*f1)
    thing.set_data(times, data1)
    return thing,

def animate(i):
    data2 = np.sin(2*np.pi*times*f2 + 4*np.pi/frames*i)
    thing.set_data(times, data2)
    return thing,


anim = animation.FuncAnimation(fig, animate, init_func = init,\
                                 frames = frames, interval = 20, blit = True)
    

anim.save("testgif.mp4", fps = 30, extra_args=['-vcodec', 'libx264'])
    

plt.show()

    