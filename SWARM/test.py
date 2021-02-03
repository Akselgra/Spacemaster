import numpy as np
import matplotlib.pyplot as plt
import aacgmv2
from SWARMprocess import SWARMprocess
from getpass import getuser
usrname = getuser()
pro = SWARMprocess()
import os
os.environ["CDF_LIB"] = "/home/" + usrname +  "/Libraries/cdf/cdf36_3-dist/lib"
from spacepy import pycdf

file = "/home/" + usrname + "/Documents/Master/Swarm_Data/SW_PREL_PCPA_PATCH_2_20131224T000000_20131224T235959_0101.cdf"

testy = pycdf.CDF(file)

print(testy["glats"].attrs)
