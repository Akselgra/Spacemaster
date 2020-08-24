"""
Requires cdf library.
"""
import numpy as np
import matplotlib.pyplot as plt
#setting path to cdf library
from getpass import getuser
usrname = getuser()
import os
os.environ["CDF_LIB"] = "/home/" + usrname +  "/Libraries/cdf/cdf/lib"
from spacepy import pycdf

cdf = pycdf.CDF("Data/testdata/testdata.cdf")
Ne = (cdf["Ne"])
latitude = np.array(cdf["Latitude"])
longitude = np.array(cdf["Longitude"])
time = (cdf["Timestamp"])

plt.plot(time, Ne)
plt.show()
