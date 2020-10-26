"""
Currently only runs on Aksels laptop.
"""
import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
import test_SWARMprocess
#setting path to cdf library
from getpass import getuser
usrname = getuser()
import os
os.environ["CDF_LIB"] = "/home/" + usrname +  "/Libraries/cdf/cdf36_3-dist/lib"
from spacepy import pycdf
from scipy.stats import pearsonr


class MovingWindow(SWARMprocess):
    """
    Class for moving a window through data and calculating cross spectrum.
    """

    def __init__(self, N = int(1e5)):
        #Initializing data
        data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"

        cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
        cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
        cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"

        self.cdfA = pycdf.CDF(cdfA_path)
        self.cdfB = pycdf.CDF(cdfB_path)
        self.cdfC = pycdf.CDF(cdfC_path)
        #Retrieving data from CDF files.
        self.NeA = self.cdfA["Ne"][:N]
        self.NeB = self.cdfB["Ne"][:N]
        self.NeC = self.cdfC["Ne"][:N]

        self.longA = self.cdfA["Longitude"][:N]
        self.longB = self.cdfB["Longitude"][:N]
        self.longC = self.cdfC["Longitude"][:N]

        self.latA = self.cdfA["Latitude"][:N]
        self.latB = self.cdfB["Latitude"][:N]
        self.latC = self.cdfC["Latitude"][:N]

        self.radA = self.cdfA["Radius"][:N]
        self.radB = self.cdfB["Radius"][:N]
        self.radC = self.cdfC["Radius"][:N]

        #Setting time to seconds after midnight
        self.seconds = self.stamp_to_sec(self.cdfA["Timestamp"][:N])
        self.stamps = self.cdfA["Timestamp"][:N]

        self.fs = 2
        
