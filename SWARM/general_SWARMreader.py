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

        self.samelength = True

        if len(self.cdfA["Ne"]) == len(self.cdfB["Ne"]) == len(self.cdfC["Ne"]):
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

            #self.fs = 2
            self.fs = 1/(self.seconds[1] - self.seconds[0])

            self.BA_shift = self.timeshift_latitude(self.latB, self.latA)
            self.BC_shift = self.timeshift_latitude(self.latB, self.latC)

        else:
            self.samelength = False
    def histmake(self, n = 100, t0 = 0, t1 = 85000, minfreq = 0, maxfreq = True,\
                 bins = 10, abs = False, norm = True):
        """
        make histograms of relative difference in integrated fouriers.

        Parameters:
            n - int; number of indices in time window used in fft_time_integral
            t0 - float; start time
            t1 - float; end time
            minfreq - float; lower integral limit
            maxfreq - float; higher integral limit
            bins - int or list; if int, number of bins. if list, bin edges.
            abs - bool; if True, relative diff will be absolute valued.
            norm - bool; if True, relative diff will be divided by highest value.

        returns:
            hists; list of histograms [BC, BA, AC]
            bins; list of bins [BC, BA, AC]
        """
        assert self.samelength, "Data sets must be without holes."
        if maxfreq:
            maxfreq = self.fs/2
        ind1 = int(self.fs*t0)
        ind2 = int(self.fs*t1)

        if len(self.NeA) < (ind2 - ind1): #fixing random, equal holes.
            ind2 = int((len(self.NeA)-ind1)*0.99)

        #timeshifting data
        NeB = self.NeB[ind1:ind2]
        NeA = self.NeA[ind1 + self.BA_shift:ind2 + self.BA_shift]
        NeC = self.NeC[ind1 + self.BC_shift:ind2 + self.BC_shift]

        #calculating integrated fourier times
        times, fft_intA = self.fft_time_integral(NeA, n, self.fs, minfreq, maxfreq)
        times, fft_intB = self.fft_time_integral(NeB, n, self.fs, minfreq, maxfreq)
        times, fft_intC = self.fft_time_integral(NeC, n, self.fs, minfreq, maxfreq)

        #finding relative difference
        fft_diff_BA = self.relative_diff(fft_intB, fft_intA, abs, norm)
        fft_diff_BC = self.relative_diff(fft_intB, fft_intC, abs, norm)
        fft_diff_AC = self.relative_diff(fft_intA, fft_intC, abs, norm)

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


    def lat_finder(self, lat1, lat0 = 0, pole = "both"):
        """
        Splits data set in high and low latitude regions
        split by lat
        returns:
            indisk - list; list on the form [[A0, A1],[B0, B1],[C0,C1]]
                           where A0...C1 are arrays of indices
                           with 1 being the region between lat0 and lat1
                           and 0 being the region between lat1 and lat0
        """

        if pole == "both":
            lowerA = np.abs(self.latA) < lat1
            higherA = np.abs(self.latA) > lat0
            is_poleA = lowerA * higherA

            lowerB = np.abs(self.latB) < lat1
            higherB = np.abs(self.latB) > lat0
            is_poleB = lowerB * higherB

            lowerC = np.abs(self.latC) < lat1
            higherC = np.abs(self.latC) > lat0
            is_poleC = lowerC * higherC

        elif pole == "north":
            lowerA = (self.latA) < lat1
            higherA = (self.latA) > lat0
            is_poleA = lowerA * higherA

            lowerB = (self.latB) < lat1
            higherB = (self.latB) > lat0
            is_poleB = lowerB * higherB

            lowerC = (self.latC) < lat1
            higherC = (self.latC) > lat0
            is_poleC = lowerC * higherC

        elif pole == "south":
            lowerA = (self.latA) > lat1
            higherA = (self.latA) < lat0
            is_poleA = lowerA * higherA

            lowerB = (self.latB) > lat1
            higherB = (self.latB) < lat0
            is_poleB = lowerB * higherB

            lowerC = (self.latC) > lat1
            higherC = (self.latC) < lat0
            is_poleC = lowerC * higherC




        # is_poleA = np.logical_not(lat0 < np.abs(self.latA) < lat)
        # is_poleB = np.logical_not(lat0 < np.abs(self.latB) < lat)
        # is_poleC = np.logical_not(lat0 < np.abs(self.latC) < lat)



        high_lat_A = np.where(is_poleA == 1)
        high_lat_B = np.where(is_poleB == 1)
        high_lat_C = np.where(is_poleC == 1)

        low_lat_A = np.where(is_poleA == 0)
        low_lat_B = np.where(is_poleB == 0)
        low_lat_C = np.where(is_poleC == 0)

        indsA = [low_lat_A, high_lat_A]
        indsB = [low_lat_B, high_lat_B]
        indsC = [low_lat_C, high_lat_C]

        indisk = [indsA, indsB, indsC]

        return indisk



    def lat_hist(self, n = 100, minfreq = 0, maxfreq = True,\
                 bins = 10, abs = False, norm = True, lat_limit = 75, lat0 = 0,\
                 pole = "both"):
        """
        makes 2 histograms for each combination of satellites
        by splitting data into high and low latitude.
        returns:
            hists - list; list of histograms. On the form [BA_high, BA_low, ... AC_low]
            bins - list; list of bins. On the form [BA_high, BA_low, ... AC_low]
        """

        assert self.samelength, "Data sets must be without holes."
        if maxfreq:
            maxfreq = self.fs/2

        #splitting by latitude
        indisk = self.lat_finder(lat_limit, lat0, pole = pole)
        high_NeA = self.NeA[indisk[0][1]]
        low_NeA = self.NeA[indisk[0][0]]

        high_NeB = self.NeB[indisk[1][1]]
        low_NeB = self.NeB[indisk[1][0]]

        high_NeC = self.NeC[indisk[2][1]]
        low_NeC = self.NeC[indisk[2][0]]

        ind1_high = 0
        ind2_high = np.min(np.array([len(high_NeA),\
                                     len(high_NeB),\
                                     len(high_NeC)]))\
                  - np.max(np.array([self.BA_shift, self.BC_shift]))

        ind1_low = 0
        ind2_low = np.min(np.array([len(low_NeA),\
                                     len(low_NeB),\
                                     len(low_NeC)]))\
                  - np.max(np.array([self.BA_shift, self.BC_shift]))

        #timeshifting data
        high_NeA = high_NeA[ind1_high + self.BA_shift:ind2_high + self.BA_shift]
        low_NeA = low_NeA[ind1_low + self.BA_shift:ind2_low + self.BA_shift]

        high_NeB = high_NeB[ind1_high:ind2_high]
        low_NeB = low_NeB[ind1_low:ind2_low]

        high_NeC = high_NeC[ind1_high + self.BC_shift:ind2_high + self.BC_shift]
        low_NeC = low_NeC[ind1_low + self.BC_shift:ind2_low + self.BC_shift]



        #calculating integrated fourier times
        times, fft_intA_high = self.fft_time_integral(high_NeA, n, self.fs, minfreq, maxfreq)
        times, fft_intA_low = self.fft_time_integral(low_NeA, n, self.fs, minfreq, maxfreq)

        times, fft_intB_high = self.fft_time_integral(high_NeB, n, self.fs, minfreq, maxfreq)
        times, fft_intB_low = self.fft_time_integral(low_NeB, n, self.fs, minfreq, maxfreq)

        times, fft_intC_high = self.fft_time_integral(high_NeC, n, self.fs, minfreq, maxfreq)
        times, fft_intC_low = self.fft_time_integral(low_NeC, n, self.fs, minfreq, maxfreq)

        #finding relative difference
        fft_diff_BA_high = self.relative_diff(fft_intB_high, fft_intA_high, abs, norm)
        fft_diff_BA_low = self.relative_diff(fft_intB_low, fft_intA_low, abs, norm)

        fft_diff_BC_high = self.relative_diff(fft_intB_high, fft_intC_high, abs, norm)
        fft_diff_BC_low = self.relative_diff(fft_intB_low, fft_intC_low, abs, norm)

        fft_diff_AC_high = self.relative_diff(fft_intA_high, fft_intC_high, abs, norm)
        fft_diff_AC_low = self.relative_diff(fft_intA_low, fft_intC_low, abs, norm)


        #making histograms
        histBA_high, binsBA_high = np.histogram(fft_diff_BA_high, bins = bins)
        histBA_low, binsBA_low = np.histogram(fft_diff_BA_low, bins = bins)

        histBC_high, binsBC_high = np.histogram(fft_diff_BC_high, bins = bins)
        histBC_low, binsBC_low = np.histogram(fft_diff_BC_low, bins = bins)

        histAC_high, binsAC_high = np.histogram(fft_diff_AC_high, bins = bins)
        histAC_low, binsAC_low = np.histogram(fft_diff_AC_low, bins = bins)


        #fixing bins
        binsBA_high = (binsBA_high[1:]+binsBA_high[:-1])/2
        binsBA_low = (binsBA_low[1:]+binsBA_low[:-1])/2

        binsBC_high = (binsBC_high[1:]+binsBC_high[:-1])/2
        binsBC_low = (binsBC_low[1:]+binsBC_low[:-1])/2

        binsAC_high = (binsAC_high[1:]+binsAC_high[:-1])/2
        binsAC_low = (binsAC_low[1:]+binsAC_low[:-1])/2


        hists = [histBA_high, histBA_low, histBC_high, histBC_low, histAC_high, histAC_low]
        bins = [binsBA_high, binsBA_low, binsBC_high, binsBC_low, binsAC_high, binsAC_low]


        return(hists, bins)







if __name__ == "__main__":
    day = "11"
    data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"
    cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501.CDF/SW_OPER_EFIA_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501_MDR_EFI_LP.cdf"
    cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501.CDF/SW_OPER_EFIB_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501_MDR_EFI_LP.cdf"
    cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501.CDF/SW_OPER_EFIC_LP_1B_201312"+day+"T000000_201312"+day+"T235959_0501_MDR_EFI_LP.cdf"
    object = GenSWARMread(cdfA_path, cdfB_path, cdfC_path)
    indisk = object.lat_finder(-90, -75, pole = "south")
    pole = object.latA[indisk[0][1]]
    times = object.seconds[indisk[0][1]]
    print(indisk[0][1])
    plt.plot(times, pole, "o")
    plt.plot(object.seconds[indisk[0][0]], object.latA[indisk[0][0]], "ro")
    plt.show()
