import numpy as np
import matplotlib.pyplot as plt

def x_1(y):
    return (0 - 2*y)
def x_2(y):
    return(2*np.pi - 2*y)


ys = np.linspace(-4*np.pi, 4*np.pi, 1000)
x1s = x_1(ys)
x2s = x_2(ys)

plt.plot(x1s, ys)
plt.plot(x2s, ys)
plt.legend(["nedre grense", "øvre grense"])
plt.grid("on")
plt.axis([-2*np.pi - 0.5, 2*np.pi + 0.5, -2*np.pi - 0.5, 2*np.pi + 0.5])
plt.xticks([-2*np.pi, -np.pi, 0, np.pi, 2*np.pi], ["$-2\pi$", "$-\pi$", "0", "$\pi$", "$2\pi$"])
plt.yticks([-2*np.pi, -np.pi, 0, np.pi, 2*np.pi], ["$-2\pi$", "$-\pi$", "0", "$\pi$", "$2\pi$"])
plt.xlabel("x")
plt.ylabel("y")
plt.title("definisjonsområde for cos(x + 2y)")
plt.show()
