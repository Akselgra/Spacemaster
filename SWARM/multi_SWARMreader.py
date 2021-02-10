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
from scipy.io import savemat
from scipy.io import loadmat


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
        self.pro = SWARMprocess()
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

    def writefile(self):
        """
        Shapes CDF files and writes them as .mat files.
        Hardcoded for december 2013.
        """
        #testing

        for i in range(self.init_loop_index, self.end_loop_index):
            files = self.gen_filenames(i)
            data = GenSWARMread(files[0], files[1], files[2], compare = True)
            # data.mlats()
            NA = len(data.NeA); NB = len(data.NeB); NC = len(data.NeC);
            N = np.min([NA, NB, NC])


            NeA = data.NeA
            NeB = data.NeB
            NeC = data.NeC

            longA = data.longA
            longB = data.longB
            longC = data.longC

            latA = data.latA
            latB = data.latB
            latC = data.latC

            radA = data.radA
            radB = data.radB
            radC = data.radC

            velA = data.velA
            velB = data.velB
            velC = data.velC

            altA = data.altA
            altB = data.altB
            altC = data.altC

            secondsA = data.secondsA
            secondsB = data.secondsB
            secondsC = data.secondsC

            stampsA = data.stampsA
            stampsB = data.stampsB
            stampsC = data.stampsC

            # mlatA = data.mlatA
            # mlatB = data.mlatB
            # mlatC = data.mlatC
            #
            # mlongA = data.mlongA
            # mlongB = data.mlongB
            # mlongC = data.mlongC
            #
            # mltA = data.mltA
            # mltB = data.mltB
            # mltC = data.mltC

            fs = data.fs

            velA = np.reshape(velA, len(velA))
            velB = np.reshape(velB, len(velB))
            velC = np.reshape(velC, len(velC))


            #equalizing holes
            NeA, NeB, NeC = self.pro.equalizer(NeA, NeB, NeC, secondsA, secondsB, secondsC, fs)
            longA, longB, longC = self.pro.equalizer(longA, longB, longC, secondsA, secondsB, secondsC, fs)
            latA, latB, latC = self.pro.equalizer(latA, latB, latC, secondsA, secondsB, secondsC, fs)
            radA, radB, radC = self.pro.equalizer(radA, radB, radC, secondsA, secondsB, secondsC, fs)
            velA, velB, velC = self.pro.equalizer(velA, velB, velC, secondsA, secondsB, secondsC, fs)
            altA, altB, altC = self.pro.equalizer(altA, altB, altC, secondsA, secondsB, secondsC, fs)
            stampsA, stampsB, stampsC = self.pro.equalizer(stampsA, stampsB, stampsC, secondsA, secondsB, secondsC, fs)
            mlatA, mlatB, mlatC = self.pro.equalizer(mlatA, mlatB, mlatC, secondsA, secondsB, secondsC, fs)
            mlongA, mlongB, mlongC = self.pro.equalizer(mlongA, mlongB, mlongC, secondsA, secondsB, secondsC, fs)
            mltA, mltB, mltC = self.pro.equalizer(mltA, mltB, mltC, secondsA, secondsB, secondsC, fs)
            secondsA, secondsB, secondsC = self.pro.equalizer(secondsA, secondsB, secondsC, secondsA, secondsB, secondsC, fs)



            day = self.day0 + i
            if day < 10:
                day = "0%g" % day
            else:
                day = "%g" % day

            path = "/home/" + usrname +  "/Documents/Master/Spacemaster/SWARM/Data/matfiles"
            file = path + "/201312" + day + ".mat"

            if not os.path.exists(path):
                os.makedirs(path)

            mdic = {"NeA":NeA, "NeB":NeB, "NeC":NeC,\
                    "longA":longA, "longB":longB, "longC":longC,\
                    "latA":latA, "latB":latB, "latC":latC,\
                    "radA":radA, "radB":radB, "radC":radC,\
                    "velA":velA, "velB":velB, "velC":velC,\
                    "altA":altA, "altB":altB, "altC":altC,\
                    #"stampsA":stampsA, "stampsB":stampsB, "stampsC":stampsC,\
                    "mlatA":mlatA, "mlatB":mlatB, "mlatC":mlatC,\
                    "mlongA":mlongA, "mlongB":mlongB, "mlongC":mlongC,\
                    "mltA":mltA, "mltB":mltB, "mltC":mltC,\
                    "secondsA":secondsA, "secondsB":secondsB, "secondsC":secondsC}
            savemat(file, mdic)




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

            # print(i)
            if data.samelength != True: #checks that data files are complete
                # print("data file %g was not complete" % i)
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
                 bins_ = 10, abs = False, norm = True, lat1 = 75, lat0 = 0,\
                 pole = "both", mlat = False):
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
            lat1 - float; high latitude limit
            lat0 - float; low latitude limit
            pole - string; "north", "south" or "both". Decides what pole we are looking at.
            mlat - bool; if True, uses magnetic latitudes

        returns:
            hists; list of histograms [BA_high, BA_low, BC_high, BC_low, AC_high, AC_low]
            bins; list of bins [BA_high, BA_low, BC_high, BC_low, AC_high, AC_low]
        """

        if norm == True:
            binlist = np.linspace(-1, 1, bins_)

        else:
            # binlist = np.linspace(-1000000, 1000000, bins_)
            binlist = np.linspace(-250000, 250000, bins_)

        for i in range(self.init_loop_index, self.end_loop_index):
            files = self.gen_filenames(i)
            data = GenSWARMread(files[0], files[1], files[2])

            #print(i)
            self.samelength = True
            if data.samelength != True: #checks that data files are complete
                #print("data file %g was not complete" % i)
                self.samelength = data.samelength
                if (self.end_loop_index - self.init_loop_index) == 1:
                    return 0, 0
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
                                        lat_limit = lat1, lat0 = lat0, pole = pole,\
                                        mlat = mlat)

            histsBA_low = histsBA_low + temp_hists[1]
            histsBA_high = histsBA_high + temp_hists[0]
            histsBC_low = histsBC_low + temp_hists[3]
            histsBC_high = histsBC_high + temp_hists[2]
            histsAC_low = histsAC_low + temp_hists[5]
            histsAC_high = histsAC_high + temp_hists[4]

            hists = np.array([histsBA_high, histsBA_low, histsBC_high, histsBC_low, histsAC_high, histsAC_low])
            self.BA_shift = data.BA_shift
            self.BC_shift = data.BC_shift
        return(hists, bins)

    def freq_sig(self, df = 0.1, n = 100,\
                 bins_ = 10, abs = False, norm = True, f0 = 0, f1 = 1):
        """
        Calculates standard deviation and mean as a function of frequency
        returns:
            freq0s - nparray; array containing lower integral limits
            sigmas - nparray; array on the form [sigmasBA, sigmasBC, sigmasAC]
            means - nparray; array on the form [meansBA, meansBC, meansAC]
        """
        N = int((f1 - f0)/df)
        freq0s = np.linspace(f0, f1-df, N)
        sigmas = np.zeros((3, N))
        means = np.zeros_like(sigmas)
        for i in range(len(freq0s)):
            print(freq0s[i])
            minfreq = freq0s[i]
            maxfreq = freq0s[i]+df
            hists, bins = self.multi_histmake(n = n, minfreq = minfreq , maxfreq = maxfreq,\
                         bins_ = bins_, abs = abs, norm = norm)
            for j in range(len(hists)):
                d_bins = bins[j][1] - bins[j][0]
                hists[j] = hists[j]/np.sum(hists[j]*d_bins)
            BAstd, BAmean = self.pro.std_mean(hists[0], bins[0])
            BCstd, BCmean = self.pro.std_mean(hists[1], bins[1])
            ACstd, ACmean = self.pro.std_mean(hists[2], bins[2])

            sigmas[0][i] = BAstd
            sigmas[1][i] = BCstd
            sigmas[2][i] = ACstd
            means[0][i] = BAmean
            means[1][i] = BCmean
            means[2][i] = ACmean


        return freq0s, sigmas, means




    def freq_sig_lat(self,df = 0.1, n = 100,\
                 bins_ = 10, abs = False, norm = True, lat1 = 75, lat0 = 0,
                 f0 = 0, f1 = 1, pole = "both", mlat = False):
        """
        Calculates standard deviation and mean as a function of frequency
        returns:
            freq0s - nparray; array containing lower integral limits
            sigmas - nparray; array on the form [sigmasBA, sigmasBC, sigmasAC]
            means - nparray; array on the form [meansBA, meansBC, meansAC]
        """
        N = int((f1 - f0)/df)
        freq0s = np.linspace(f0, f1-df, N)
        sigmas = np.zeros((6, N))
        means = np.zeros_like(sigmas)
        for i in range(len(freq0s)):
            print(freq0s[i])
            minfreq = freq0s[i]
            maxfreq = freq0s[i]+df
            hists, bins = self.multi_histmake_lat(n = n, minfreq = minfreq , maxfreq = maxfreq,\
                         bins_ = bins_, abs = abs, norm = norm, lat1 = lat1, lat0 = lat0, pole = pole, mlat = mlat)
            for j in range(len(hists)):
                d_bins = bins[j][1] - bins[j][0]
                hists[j] = hists[j]/np.sum(hists[j]*d_bins)


            BAstd_high, BAmean_high = self.pro.std_mean(hists[0], bins[0])
            BAstd_low, BAmean_low = self.pro.std_mean(hists[1], bins[1])
            BCstd_high, BCmean_high = self.pro.std_mean(hists[2], bins[2])
            BCstd_low, BCmean_low = self.pro.std_mean(hists[3], bins[3])
            ACstd_high, ACmean_high = self.pro.std_mean(hists[4], bins[4])
            ACstd_low, ACmean_low = self.pro.std_mean(hists[5], bins[5])

            sigmas[0][i] = BAstd_high
            sigmas[1][i] = BAstd_low
            sigmas[2][i] = BCstd_high
            sigmas[3][i] = BCstd_low
            sigmas[4][i] = ACstd_high
            sigmas[5][i] = ACstd_low
            means[0][i] = BAmean_high
            means[1][i] = BAmean_low
            means[2][i] = BCmean_high
            means[3][i] = BCmean_low
            means[4][i] = ACmean_high
            means[5][i] = ACmean_low


        return freq0s, sigmas, means






if __name__ == "__main__":
    object = MultiSWARM(2013, 12, 9, 2013, 12, 9)
    object.writefile()
