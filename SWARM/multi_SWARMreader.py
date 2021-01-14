"""
CDF_LIB must be in the correct file path to run
"""
import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
from general_SWARMreader import GenSWARMread
#setting path to cdf library
from getpass import getuser
usrname = getuser()
import os
os.environ["CDF_LIB"] = "/home/" + usrname +  "/Libraries/cdf/cdf36_3-dist/lib"
from spacepy import pycdf


class MultiSWARM():
    """
    Class for reading multiple files of SWARM data.
    Designed for specific file structures.
    """

    def __init__(self, year1, month1, day1, year2, month2, day2,\
                 day0 = 9, month0 = 12, year0 = 2013 ):
        """
        parameters:
            year1 - int; year of first data file
            month1 - int; month of first data file
            day1 - int; day of first data file
            year2 - int; year of last data file
            month2 - int; month of last data file
            day2 - int; day of last data file
            day0 - int; day of first file
            month0 - int; month of first file
            year0 - int; year of first file

        """
        self.year0 = year0
        self.month0 = month0
        self.day0 = day0
        self.year1 = year1
        self.month1 = month1
        self.day1 = day1
        self.year2 = year2
        self.month2 = month2
        self.day2 = day2

        self.data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"

        #calculating indices for looping over data files
        start_days = 0
        days_per_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        tot_month_days0 = 0
        tot_month_days1 = 0
        tot_month_days2 = 0
        for i in range(month0-1):
            tot_month_days0 += days_per_month[i]

        for i in range(month1-1):
            tot_month_days1 += days_per_month[i]

        for i in range(month2-1):
            tot_month_days2 += days_per_month[i]

        init_loop_index = 365*year1 + tot_month_days1 + day1\
                        -(365*year0 + tot_month_days0 + day0)

        end_loop_index = 365*year2 + tot_month_days2 + day2\
                       - (365*year0 + tot_month_days0 + day0)

        self.init_loop_index = init_loop_index
        self.end_loop_index = end_loop_index + 1

        if self.init_loop_index < 0:
            raise ValueError("initial date must be higher or equal to base date")
        if self.end_loop_index < self.init_loop_index:
            raise ValueError("end date must be higher than initial date")






    def filenames(self, year, month, day):
        """
        Takes a date and returns a list of file path strings.
        parameters:
            year - int; year of data files yyyy
            month - int; month of data files mm
            day - int; day of data files dd
        returns:
            files - list; file paths on the form ["satA", "satB", "satC"]
        """
        data_path = self.data_path

        year = str(year)
        if len(year) != 4:
            raise ValueError("year must be on form yyyy")

        month = str(month)
        if len(month) == 1:
            month = "0" + month

        day = str(day)
        if len(day) == 1:
            day = "0" + day

        date = year + month + day
        t0 = "T000000_"
        t1 = "T235959_0501"

        cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_" + date + t0 + date\
        + t1 + ".CDF/SW_OPER_EFIA_LP_1B_" + date + t0 + date + t1 + "_MDR_EFI_LP.cdf"

        cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_" + date + t0 + date\
         + t1 + ".CDF/SW_OPER_EFIB_LP_1B_" + date + t0 + date + t1 + "_MDR_EFI_LP.cdf"

        cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_" + date + t0 + date\
         + t1 + ".CDF/SW_OPER_EFIC_LP_1B_" + date + t0 + date + t1 + "_MDR_EFI_LP.cdf"

        files = [cdfA_path, cdfB_path, cdfC_path]

        return(files)

    def gen_filenames(self, ind):
        """
        Generalized version of filenames. Takes an index n and returns
        filepaths of sats A, B and C for the day of the given index in directory.
        example: ind 0 gives first day. ind 21 gives 21st day.
        parameters:
            ind - int; index of day in directory
        returns:
            files - list; list of filepaths on the form [satA, satB, satC]
        """

        pathA = self.data_path + "/Sat_A"
        dirsA = os.listdir(pathA)
        dirsA.sort()
        pathA = pathA + "/" + dirsA[ind]
        cdfA = os.listdir(pathA)
        cdfA.sort()
        pathA = pathA + "/" + cdfA[1]

        pathB = self.data_path + "/Sat_B"
        dirsB = os.listdir(pathB)
        dirsB.sort()
        pathB = pathB + "/" + dirsB[ind]
        cdfB = os.listdir(pathB)
        cdfB.sort()
        pathB = pathB + "/" + cdfB[1]

        pathC = self.data_path + "/Sat_C"
        dirsC = os.listdir(pathC)
        dirsC.sort()
        pathC = pathC + "/" + dirsC[ind]
        cdfC = os.listdir(pathC)
        cdfC.sort()
        pathC = pathC + "/" + cdfC[1]

        files = [pathA, pathB, pathC]

        return(files)


    def multi_histmake(self, n = 100, minfreq = 0, maxfreq = True,\
                 bins_ = 10, abs = False, norm = True):
        """
        make histograms of relative difference in integrated fouriers using
        multiple dates of data

        Parameters:
            n - int; number of indices in time window used in fft_time_integral
            t0 - float; start time
            t1 - float; end time
            minfreq - float; lower integral limit
            maxfreq - float; higher integral limit
            bins_ - int; number of bin edges.

        returns:
            hists; list of histograms [BC, BA, AC]
            bins; list of bins [BC, BA, AC]
        """

        binlist = np.linspace(-1, 1, bins_)

        for i in range(self.init_loop_index, self.end_loop_index):
            files = self.gen_filenames(i)
            data = GenSWARMread(files[0], files[1], files[2])

            print(i)
            if data.samelength != True: #checks that data files are complete
                print("data file %g was not complete" % i)
                continue

            if maxfreq:
                maxfreq = data.fs/2
            histsA = np.zeros_like(binlist[:-1])
            histsB = np.zeros_like(binlist[:-1])
            histsC = np.zeros_like(binlist[:-1])
            t0 = data.seconds[200]
            t1 = data.seconds[-1600] #data set is sliced to make room for index shifting


            temp_hists, bins = data.histmake(n = n, minfreq = minfreq, maxfreq = maxfreq,\
                                        t0 = t0, t1 = t1, bins = binlist, abs = abs, norm = norm)

            histsA = histsA + temp_hists[0]
            histsB = histsB + temp_hists[1]
            histsC = histsC + temp_hists[2]

            hists = np.array([histsA, histsB, histsC])
        return(hists, bins)


    def multi_histmake_lat(self, n = 100, minfreq = 0, maxfreq = True,\
                 bins_ = 10, abs = False, norm = True, lat1 = 75, lat0 = 0):
        """
        make histograms of relative difference in integrated fouriers using
        multiple dates of data

        Parameters:
            n - int; number of indices in time window used in fft_time_integral
            t0 - flofor i in range(object.init_loop_index, object.end_loop_index):
        files = object.gen_filenames(i)
        genread = GenSWARMread(files[0], files[1], files[2])
        if genread.samelength != True:
            print("%g failed" % i)
            continue
        # binlist = np.linspace(-1, 1, 50)
        # hists, bins = genread.lat_hist(lat_limit = 85, bins = binlist)
        #
        # print(hists[1])
        # plt.figure(1)
        # plt.bar(bins[0],hists[0], width = 0.9*(bins[0][1] - bins[0][0]))
        # plt.figure(2)
        # plt.bar(bins[1], hists[1], width = 0.9*(bins[1][1] - bins[1][0]))
        # plt.show()

        indisk = genread.lat_finder(85, 80)
        plt.plot(genread.seconds[indisk[0][1]], genread.latA[indisk[0][1]], ".")
        plt.show()at; start time
            t1 - float; end time
            minfreq - float; lower integral limit
            maxfreq - float; higher integral limit
            bins_ - int; number of bin edges.
            lat1 - float; high latitude limit
            lat0 - float; low latitude limit

        returns:
            hists; list of histograms [BC, BA, AC]
            bins; list of bins [BC, BA, AC]
        """

        binlist = np.linspace(-1, 1, bins_)

        for i in range(self.init_loop_index, self.end_loop_index):
            files = self.gen_filenames(i)
            data = GenSWARMread(files[0], files[1], files[2])

            print(i)
            if data.samelength != True: #checks that data files are complete
                print("data file %g was not complete" % i)
                continue

            if maxfreq:
                maxfreq = data.fs/2
            histsBA_low = np.zeros_like(binlist[:-1])
            histsBA_high = np.zeros_like(binlist[:-1])
            histsBC_low = np.zeros_like(binlist[:-1])
            histsBC_high = np.zeros_like(binlist[:-1])
            histsAC_low = np.zeros_like(binlist[:-1])
            histsAC_high = np.zeros_like(binlist[:-1])
            t0 = data.seconds[200]
            t1 = data.seconds[-1600] #data set is sliced to make room for index shifting


            temp_hists, bins = data.lat_hist(n = n, minfreq = minfreq, maxfreq = maxfreq,\
                                        bins = binlist, abs = abs, norm = norm,\
                                        lat_limit = lat1, lat0 = lat0)

            histsBA_low = histsBA_low + temp_hists[1]
            histsBA_high = histsBA_high + temp_hists[0]
            histsBC_low = histsBC_low + temp_hists[3]
            histsBC_high = histsBC_high + temp_hists[2]
            histsAC_low = histsAC_low + temp_hists[5]
            histsAC_high = histsAC_high + temp_hists[4]

            hists = np.array([histsBA_high, histsBA_low, histsBC_high, histsBC_low, histsAC_high, histsAC_low])
        return(hists, bins)

    def testymctestface(self):
        """
        tests testy testfaces
        """
        for i in range(self.init_loop_index, self.end_loop_index):
            files = self.gen_filenames(i)
            data = GenSWARMread(files[0], files[1], files[2])
            if data.samelength != True:
                continue

            shifts = 7000
            start = 0
            stop = 20000
            lat1 = data.latA
            lat2 = data.latC
            meandist = np.zeros(shifts)
            indices = np.arange(shifts)

            for j in range(shifts):
                meandist[j] = np.mean(np.abs(lat1[start:stop] - lat2[start + j: stop + j]))

            plt.plot(indices/2,1 - meandist/np.max(meandist))
            plt.title("day %g" % (self.day0 + i))
            plt.xlabel("Time shifted [s]")
            plt.ylabel("normalized mean distance in latitude")
            plt.show()






if __name__ == "__main__":
    object = MultiSWARM(2013, 12, 9, 2013, 12, 31)
    hists, bins = object.multi_histmake_lat(minfreq = 1/3, maxfreq = 2/3, bins_ = 50, lat0 = 75, lat1 = 90)
    width = bins[0][1] - bins[0][0]
    plt.figure(1)
    plt.bar(bins[0], hists[0], width = 0.9*width)
    plt.title("high lat")
    plt.figure(2)
    plt.bar(bins[1], hists[1], width = 0.9*width)
    plt.title("low lat")
    plt.show()
    # minfreq = 0
    # maxfreq = 1/3
    # hists, bins = object.multi_histmake(bins_ = 200, minfreq = minfreq, maxfreq = maxfreq)
    #
    #
    #
    # widthA = bins[0][1] - bins[0][0]
    # BA_hist = hists[0]/np.sum(hists[0]*widthA)
    #
    # std = np.std(BA_hist*bins[0])
    # mean = np.mean(BA_hist*bins[0])
    # x = np.linspace(-1, 1, 1000)
    # pro = SWARMprocess()
    # gaussian = pro.gauss_curve(x, std = std, mean = mean)
    #
    # print(std)
    # print(mean)
    #
    # plt.figure(1)
    # plt.bar(bins[0], BA_hist, width = 0.9*widthA)
    # plt.xlabel("relative difference B - A")
    # plt.ylabel("Normalized occurence")
    # plt.title("B - A. integral limits: %g to %g" % (minfreq, maxfreq))
    # plt.plot(x, gaussian, "r")
    # plt.legend(["gaussian approximation", "B-A histogram"])
    # plt.show()
