import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from getpass import getuser
usrname = getuser()

path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"
file = path + "/test.mat"

infile = loadmat(file)
NeA = infile["NeA"][0]
secondsA = infile["secondsA"][0]
plt.plot(secondsA, NeA)
plt.show()
