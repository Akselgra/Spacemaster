"""
Currently only runs on Aksels laptop.
"""
import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
#setting path to cdf library
from getpass import getuser
usrname = getuser()
import os
os.environ["CDF_LIB"] = "/home/" + usrname +  "/Libraries/cdf/cdf36_3-dist/lib"
from spacepy import pycdf
from scipy.stats import pearsonr




class WaveInspect(SWARMprocess):
    """
    Class for inspecting waves
    """

    def __init__(self, N = int(1e5)):
        data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"
        # cdfA = pycdf.CDF("Data/Sat_A/SW_OPER_EFIA_LP_1B_20131202T101113_20131202T140109_0501.CDF/SW_OPER_EFIA_LP_1B_20131202T101113_20131202T140109_0501_MDR_EFI_LP.cdf")
        # cdfB = pycdf.CDF("Data/Sswarm-diss.eo.esa.intat_B/SW_OPER_EFIB_LP_1B_20131202T114445_20131202T153609_0501.CDF/SW_OPER_EFIB_LP_1B_20131202T114445_20131202T153609_0501_MDR_EFI_LP.cdf")
        # cdfC = pycdf.CDF("Data/Sat_C/SW_OPER_EFIC_LP_1B_20131204T094004_20131204T223759_0501.CDF/SW_OPER_EFIC_LP_1B_20131204T094004_20131204T223759_0501_MDR_EFI_LP.cdf")

        cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
        cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
        cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"

        cdfA = pycdf.CDF(cdfA_path)
        cdfB = pycdf.CDF(cdfB_path)
        cdfC = pycdf.CDF(cdfC_path)

        self.NeA = cdfA["Ne"][:N]
        self.NeB = cdfB["Ne"][:N]
        self.NeC = cdfC["Ne"][:N]

        self.longA = cdfA["Longitude"][:N]
        self.longB = cdfB["Longitude"][:N]
        self.longC = cdfC["Longitude"][:N]

        self.latA = cdfA["Latitude"][:N]
        self.latB = cdfB["Latitude"][:N]
        self.latC = cdfC["Latitude"][:N]

        self.radA = cdfA["Radius"][:N]
        self.radB = cdfB["Radius"][:N]
        self.radC = cdfC["Radius"][:N]

        self.seconds = self.stamp_to_sec(cdfA["Timestamp"][:N])

        self.pole_finder()


    def pole_finder(self):
        """
        finds the indices where the satellites
        are between 80 and 90 degrees latitude
        """
        init_cut = 5000
        is_poleA = np.logical_not(self.latA[:init_cut] < 80)
        is_poleB = np.logical_not(self.latB[:init_cut] < 80)
        is_poleC = np.logical_not(self.latC[:init_cut] < 80)

        self.indA = np.where(is_poleA == 1)
        self.indB = np.where(is_poleB == 1)
        self.indC = np.where(is_poleC == 1)



if __name__ == "__main__":
    object = WaveInspect()
    object.pole_finder()
