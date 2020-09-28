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
        #Initializing data
        data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"

        cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
        cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
        cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"

        cdfA = pycdf.CDF(cdfA_path)
        cdfB = pycdf.CDF(cdfB_path)
        cdfC = pycdf.CDF(cdfC_path)
        #Retrieving data from CDF files.
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

        #Setting time to seconds after midnight
        self.seconds = self.stamp_to_sec(cdfA["Timestamp"][:N])
        self.stamps = cdfA["Timestamp"][:N]
        #finding first indices of pole region
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


    def pole_case_study(self):
        """
        Performs a case study of electron density movement
        in the polar region.
        """

        NeA = self.NeA[self.indB]
        NeB = self.NeB[self.indB]
        NeC = self.NeC[self.indB]
        seconds = self.seconds[self.indB]

        plt.plot(seconds, NeB)
        plt.plot(seconds, NeA)
        plt.plot(seconds, NeC)
        plt.xlabel("Seconds since midnight")
        plt.ylabel("Electron density")
        plt.title("Electron density measurements at polar region")
        plt.legend(["Sat B", "Sat A", "Sat C"])
        plt.show()

        AB_pearson = pearsonr(NeB, NeA)[0]
        print("The pearson correlation coefficient between A and B is %g" % \
                                                                    AB_pearson)
        start = self.indB[0][0]
        stop = self.indB[0][-1]

        BA_shift = self.timeshift_latitude(self.latB, self.latA, start,\
         stop, shifts = 7500)

        BC_shift = self.timeshift_latitude(self.latB, self.latC, start,\
         stop, shifts = 7500)

        NeA = self.NeA[start + BA_shift:stop + BA_shift +1]
        NeC = self.NeC[start + BC_shift:stop + BC_shift +1]

        plt.plot(seconds, NeB)
        plt.plot(seconds, NeA)
        plt.plot(seconds, NeC)
        plt.xlabel("Seconds since midnight")
        plt.ylabel("Electron density [cm⁻¹]")
        plt.title("Shifted electron density measurements")
        plt.legend(["Sat B", "Sat A", "Sat C"])
        plt.show()

        AB_shift_pearson = pearsonr(NeB, NeA)[0]
        print("After shifting, the pearson correlation coefficient between A\
         and B is %g" % AB_shift_pearson)





if __name__ == "__main__":
    object = WaveInspect()
    object.pole_case_study()
