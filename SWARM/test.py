import numpy as np
import matplotlib.pyplot as plt


thing = np.random.randint(0, 2, size = 10)
print(thing)

liste = []

for i in thing:
    if i == 0:
        liste.append("k")
    else:
        liste.append("m")

print(np.sum(thing))