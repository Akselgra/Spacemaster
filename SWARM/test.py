import numpy as np

x = np.array([True, True, False])
y = np.array([False, False, False])

f = np.array([1,2,3])
x = 1 < f
y = f < 3
print(x)
print(y)
print(x*y)
