"""
CDF_LIB must be in the correct file path to run
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


class GenSWARMread(SWARMprocess):
    """
    Generalized class for reading SWARM data
    """

    def __init__(self, filenameA, filenameB, filenameC, N = True):
        """
        Constructor
        parameters:
            filenameA - str; path to and name of CSD A file
            filenameB - str; path to and name of CSD B file
            filenameC - str; path to and name of CSD C file
        """

        self.cdfA = pycdf.CDF(filenameA)
        self.cdfB = pycdf.CDF(filenameB)
        self.cdfC = pycdf.CDF(filenameC)
        #Retrieving data from CDF files.
        if N:
            N = len(self.cdfA["Ne"])
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

    def histmake(self, n = 100, t0 = 0, t1 = 85000, minfreq = 0, maxfreq = True,\
                 bins = 10):
        """
        make histograms of relative difference in integrated fouriers.

        Parameters:
            n - int; number of indices in time window used in fft_time_integral
            t0 - float; start time
            t1 - float; end time
            minfreq - float; lower integral limit
            maxfreq - float; higher integral limit

        returns:
            hists; list of histograms [BC, BA, AC]
            bins; list of bins [BC, BA, AC]
        """
        if maxfreq:
            maxfreq = self.fs/2
        ind1 = int(self.fs*t0)
        ind2 = int(self.fs*t1)

        #timeshifting data
        NeB = self.NeB[ind1:ind2]
        NeA = self.NeA[ind1 + self.BA_shift:ind2 + self.BA_shift]
        NeC = self.NeC[ind1 + self.BC_shift:ind2 + self.BC_shift]

        #calculating integrated fourier times
        times, fft_intA = self.fft_time_integral(NeA, n, self.fs, minfreq, maxfreq)
        times, fft_intB = self.fft_time_integral(NeB, n, self.fs, minfreq, maxfreq)
        times, fft_intC = self.fft_time_integral(NeC, n, self.fs, minfreq, maxfreq)

        #finding relative difference
        fft_diff_BA = self.relative_diff(fft_intB, fft_intA)
        fft_diff_BC = self.relative_diff(fft_intB, fft_intC)
        fft_diff_AC = self.relative_diff(fft_intB, fft_intC)

        #making histograms
        histBA, binsBA = np.histogram(fft_diff_BA, bins = bins)
        histBC, binsBC = np.histogram(fft_diff_BC, bins = bins)
        histAC, binsAC = np.histogram(fft_diff_AC, bins = bins)

        #fixing bins
        binsBA = (binsBA[1:]+binsBA[:-1])/2
        binsBC = (binsBC[1:]+binsBC[:-1])/2
        binsAC = (binsAC[1:]+binsAC[:-1])/2

        hists = [histBA, histBC, histAC]
        bins = [binsBA, binsBC, binsAC]

        return(hists, bins)





if __name__ == "__main__":
    day = "29"
    data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"
    cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501.CDF/SW_OPER_EFIA_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501_MDR_EFI_LP.cdf"
    cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501.CDF/SW_OPER_EFIB_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501_MDR_EFI_LP.cdf"
    cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501.CDF/SW_OPER_EFIC_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501_MDR_EFI_LP.cdf"
    object = GenSWARMread(cdfA_path, cdfB_path, cdfC_path)
    hists, bins = object.histmake(minfreq = 0, maxfreq = 0.33)
    BA_hist1 = hists[0]
    BA_bins1 = bins[0]
    width1 = BA_bins1[1] - BA_bins1[0]

    hists, bins = object.histmake(minfreq = 0.33, maxfreq = 0.66)
    BA_hist2 = hists[0]
    BA_bins2 = bins[0]
    width2 = BA_bins2[1] - BA_bins2[0]

    hists, bins = object.histmake(minfreq = 0.66, maxfreq = 1)
    BA_hist3 = hists[0]
    BA_bins3 = bins[0]
    width3 = BA_bins3[1] - BA_bins3[0]

    df1 = BA_bins1[1] - BA_bins1[0]
    norm1 = np.sum(BA_hist1)*df1

    df2 = BA_bins2[1] - BA_bins2[0]
    norm2 = np.sum(BA_hist2)*df2

    df3 = BA_bins3[1] - BA_bins3[0]
    norm3 = np.sum(BA_hist3)*df3

    plt.figure(1)
    plt.bar(BA_bins1, BA_hist1/norm1, width = width1*0.9)
    plt.title("F: 0 - 0.33")
    plt.xlabel("relative difference")
    plt.ylabel("Normalized occurance rate")
    plt.savefig("Figures/testhist_F_0_033.png")

    plt.figure(2)
    plt.bar(BA_bins2, BA_hist2/norm3, width = width2*0.9)
    plt.title("F: 0.33 - 0.66")
    plt.xlabel("relative difference")
    plt.ylabel("Normalized occurence rate")
    plt.savefig("Figures/testhist_F_033_066.png")

    plt.figure(3)
    plt.bar(BA_bins3, BA_hist3/norm3, width = width3*0.9)
    plt.title("F: 0.66 - 1")
    plt.xlabel("relative difference")
    plt.ylabel("Normalized occurence rate")
    plt.savefig("Figures/testhist_F_066_1.png")
    plt.show()

    plt.plot(BA_bins1, BA_hist1/norm1)
    plt.plot(BA_bins2, BA_hist2/norm2)
    plt.plot(BA_bins3, BA_hist3/norm3)
    plt.legend(["0 - 0.33", "0.33 - 0.66", "0.66 - 1"])
    plt.show()

    BC_hist = hists[1]
    BC_bins = bins[1]
    dw = BC_bins[1] - BC_bins[0]
    norm = np.sum(BC_hist)*dw

    plt.plot(BA_bins3, BA_hist3/norm3)
    plt.plot(BC_bins, BC_hist/norm)
    plt.show()
