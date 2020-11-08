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
        Calculates cross spectral densities and fourier transforms over time.
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
        AC_CSDs = []
        B_PSDs = []
        A_PSDs = []
        C_PSDs = []


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
            AC_cross_spec = self.cross_spectrum(dataA, dataC)[:n_freq]
            A_fft = np.fft.fft(dataA)[:n_freq]/(2*n_freq)
            B_fft = np.fft.fft(dataB)[:n_freq]/(2*n_freq)
            C_fft = np.fft.fft(dataC)[:n_freq]/(2*n_freq)
            BA_CSDs.append(np.abs(BA_cross_spec))
            BC_CSDs.append(np.abs(BC_cross_spec))
            AC_CSDs.append(np.abs(AC_cross_spec))
            A_PSDs.append(np.abs(A_fft))
            B_PSDs.append(np.abs(B_fft))
            C_PSDs.append(np.abs(C_fft))

        BA_CSDs = np.array(BA_CSDs)
        BC_CSDs = np.array(BC_CSDs)
        AC_CSDs = np.array(AC_CSDs)
        A_PSDs = np.array(A_PSDs)
        B_PSDs = np.array(B_PSDs)
        C_PSDs = np.array(C_PSDs)
        self.freqs = freqs
        self.times = times
        self.BA_CSDs = BA_CSDs
        self.BC_CSDs = BC_CSDs
        self.AC_CSDs = AC_CSDs
        self.A_PSDs = A_PSDs
        self.B_PSDs = B_PSDs
        self.C_PSDs = C_PSDs
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
        plt.title("logarithmic CSD BA")
        plt.show()

        loggy2 = np.log10(self.BC_CSDs)
        plt.contourf(Times, Freqs, loggy2, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("logarithmic CSD BC")
        plt.show()

        loggy3 = np.log10(self.AC_CSDs)
        plt.contourf(Times, Freqs, loggy3, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("logarithmic CSD AC")
        plt.show()

        logA = np.log10(self.A_PSDs)
        logB = np.log10(self.B_PSDs)
        logC = np.log10(self.C_PSDs)



    def diffplot(self):
        """
        Plots comparisons between CSDs and PSDs
        """

        assert self.solved, "solver method has not been called."

        Freqs, Times = np.meshgrid(self.freqs, self.times[:-1])

        #normalizing CSDs and PSDs
        # BA_CSDs = self.BA_CSDs/np.max(self.BA_CSDs)
        # BC_CSDs = self.BC_CSDs/np.max(self.BC_CSDs)
        # AC_CSDs = self.AC_CSDs/np.max(self.AC_CSDs)
        # A_PSDs = self.A_PSDs/np.max(self.A_PSDs)
        # B_PSDs = self.B_PSDs/np.max(self.B_PSDs)
        # C_PSDs = self.C_PSDs/np.max(self.C_PSDs)

        BA_CSDs = self.BA_CSDs
        BC_CSDs = self.BC_CSDs
        AC_CSDs = self.AC_CSDs
        A_PSDs = self.A_PSDs
        B_PSDs = self.B_PSDs
        C_PSDs = self.C_PSDs

        logBA = np.log10(BA_CSDs)
        logBC = np.log10(BC_CSDs)
        logAC = np.log10(AC_CSDs)
        logA = np.log10(A_PSDs)
        logB = np.log10(B_PSDs)
        logC = np.log10(C_PSDs)


        savepath = "contour_figs/"
        timestr = "_t0 = %g_t1=%g.png" % (self.t0, self.t1)
        timetitle = ",t = [%g, %g]" % (self.t0, self.t1)
        plt.contourf(Times, Freqs, logBA, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Logarithmic CSD BA" + timetitle )
        plt.savefig(savepath + "CSD_BA" + timestr)
        plt.close()

        plt.contourf(Times, Freqs, logBC, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Logarithmic CSD BC" + timetitle )
        plt.savefig(savepath + "CSD_BC" + timestr)
        plt.close()

        plt.contourf(Times, Freqs, logAC, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Logarithmic CSD AC" + timetitle )
        plt.savefig(savepath + "CSD_AC" + timestr)
        plt.close()


        plt.contourf(Times, Freqs, logA, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Logarithmic PSD A" + timetitle)
        plt.savefig(savepath + "PSD_A" + timestr)
        plt.close()

        plt.contourf(Times, Freqs, logB, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Logarithmic PSD B" + timetitle)
        plt.savefig(savepath + "PSD_B" + timestr)
        plt.close()

        plt.contourf(Times, Freqs, logC, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Logarithmic PSD C" + timetitle)
        plt.savefig(savepath + "PSD_C" + timestr)
        plt.close()

        BA_AC_diff = logBA - logAC
        BA_BC_diff = logBA - logBC
        BC_AC_diff = logBC - logAC

        plt.contourf(Times, Freqs, BA_AC_diff, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Log10(CSD BA) - log10(CSD AC)" + timetitle)
        plt.savefig(savepath + "BA_AC_diff" + timestr)
        plt.close()

        plt.contourf(Times, Freqs, BA_BC_diff, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Log10(CSD BA) - log10(CSD BC)" + timetitle)
        plt.savefig(savepath + "BA_BC_diff" + timestr)
        plt.close()

        plt.contourf(Times, Freqs, BC_AC_diff, cmap = "magma")
        plt.xlabel("Seconds after midnight")
        plt.ylabel("Frequency")
        plt.colorbar()
        plt.title("Log10(CSD BC) - log10(CSD AC)" + timetitle)
        plt.savefig(savepath + "BC_AC_diff" + timestr)
        plt.close()


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

    t0s = [0, 1000, 1100]
    t1s = [49000, 8000, 1500]
    ns = [10000, 1000, 50]
    for i in range(len(t0s)):
        n_window = int((t1s[i]-t0s[i])/ns[i]*2)*2
        window = windows.general_gaussian(n_window, 1, sig = n_window/4)
        object.solver(t0s[i], t1s[i], ns[i], window = window)
        object.diffplot()
