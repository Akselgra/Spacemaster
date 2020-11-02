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
import scipy.signal.windows as windows


class MovingWindow(SWARMprocess):
    """
    Class for moving a window through data and calculating cross spectrum.
    """

    def __init__(self, N = int(1e5)):
        self.N = N
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

        self.BA_shift = self.timeshift_latitude(self.latB, self.latA)
        self.BC_shift = self.timeshift_latitude(self.latB, self.latC)

        self.solved = False




    def solver(self, t0, t1, n, window = 1):
        """
        Calculates cross spectral densities over time.
        __________________________________________________
        Arguments:
        t0 - initial time in seconds
        t1 - end time in seconds
        n -  number of time points
        window - window to be applied.
        """
        if int(self.fs*t1) > self.N:
            raise ValueError("end time is larger than data set")

        n += 1

        if t0 < t1/n:
            t0 = t1/n+1

        times = np.linspace(t0, t1, n)
        indices = np.zeros_like(times) #indices of windows
        for i in range(len(indices)):
            indices[i] = int(self.fs*times[i])

        n_freq = int((indices[1] - indices[0])) #length of frequency domain
        freqs = np.linspace(0, self.fs/2, n_freq)
        BA_CSDs = []
        BC_CSDs = []


        for i in range(len(indices)-1):
            ind0 = int(indices[i] - n_freq)
            ind1 = int(indices[i+1])

            n_window = ind1 - ind0
            window = windows.general_gaussian(n_window, 1, sig = n_window/8)

            dataA = self.NeA[int(ind0 + self.BA_shift):\
                             int(ind1 + self.BA_shift)]
            dataB = self.NeB[int(ind0):int(ind1)]
            dataC = self.NeC[int(ind0 + self.BC_shift):\
                             int(ind1 + self.BC_shift)]

            dataA = dataA - np.mean(dataA)
            dataB = dataB - np.mean(dataB)
            dataC = dataC - np.mean(dataC)

            dataA = dataA*window
            dataB = dataB*window
            dataC = dataC*window

            BA_cross_spec = self.cross_spectrum(dataB, dataA)[:n_freq]
            BC_cross_spec = self.cross_spectrum(dataB, dataC)[:n_freq]
            BA_CSDs.append(np.abs(BA_cross_spec))
            BC_CSDs.append(np.abs(BC_cross_spec))

        BA_CSDs = np.array(BA_CSDs)
        BC_CSDs = np.array(BC_CSDs)
        self.freqs = freqs
        self.times = times
        self.BA_CSDs = BA_CSDs
        self.BC_CSDs = BC_CSDs
        self.solved = True
        self.t0 = t0
        self.t1 = t1


    def contourplot(self):
        """
        If solve has been called, plots the results.
        """

        assert self.solved, "solver method has not been called."

        Freqs, Times = np.meshgrid(self.freqs, self.times[:-1])

        loggy = np.log10(self.BA_CSDs)
        plt.contourf(Times, Freqs, loggy, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("logarithmic CSD")
        plt.show()

        plt.plot(self.seconds[int(2*self.t0):int(2*self.t1)],\
         self.NeB[int(2*self.t0):int(2*self.t1)])
        plt.show()
        # plt.plot(self.seconds, self.NeB)
        # plt.show()
if __name__ == "__main__":
    object = MovingWindow()
    # t0 = 3750
    # t1 = 4600
    # n = 100

    # t0 = 3000
    # t1 = 6000
    # n = 1000

    # t0 = 0
    # t1 = 49000
    # n = 10000

    t0s = [0, 3000, 3750]
    t1s = [49000, 8000, 4600]
    ns = [10000, 1000, 100]
    for i in range(len(t0s)):
        n_window = int((t1s[i]-t0s[i])/ns[i]*2)*2
        window = windows.general_gaussian(n_window, 1, sig = n_window/4)
        object.solver(t0s[i], t1s[i], ns[i], window = window)
        object.contourplot()
