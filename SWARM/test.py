from scipy.io import loadmat
from scipy.io import savemat
import numpy as np


infile = loadmat("Data/matfiles/20131209.mat")
print(infile["NeA"])