import numpy as np
import matplotlib.pyplot as plt

def f(y):
    return (1 - y**2)

ys = np.linspace(-100, 100, 1000)
xs = f(ys)

plt.plot(xs, ys)
plt.xlabel("x")
plt.ylabel("y")
plt.show()