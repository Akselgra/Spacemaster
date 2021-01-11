import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from getpass import getuser
usrname = getuser()

path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"
dirsA = listdir(path + "/Sat_A")
cdfA = listdir(path + "/Sat_A/" + dirsA[53])
cdfA_path = path + "/sat_A/" + dirsA[53] + "/" + cdfA[0]
print(cdfA_path)
