from SWARMprocess import SWARMprocess
import numpy as np
import matplotlib.pyplot as plt
from general_matreader import MatReader
import os
from scipy.io import loadmat
from scipy.io import savemat



class MultiMat():
    """
    class for reading multiple .mat files of SWARM data.
    Designed for December 2013
    """

    def __init__(self, day0, day1, day_0 = 9, day_1 = 31,\
                 data_path = "Data/matfiles"):
        """
        Parameters:
                day0 - initial day. Must be between day_0 and day_1.
                day1 - final day. Must be between day_0 and day_1.
                day_0 - first day in dataset
                day_1 - last day in dataset
        """
        assert (day0 >= day_0 and day0 <= day_1), "day0 must be between %g and %g." % (day_0, day_1)
        assert (day1 >= day_0 and day1 <= day_1), "day1 must be between %g and %g" % (day_0, day_1)
        assert (day1 >= day0), "day1 must be larger or equal to day0"

        self.day_0 = day_0

        self.data_path = data_path #path to .mat files

        self.init_loop_index = day0 - day_0
        self.end_loop_index = day1 - day_0 + 1
        self.pro = SWARMprocess()


    def filepath(self, ind):
        """
        Takes an index ind and returns filepath for day of the given index in
        directory.
        Example: if ind is 9, returns filename of 9th of December.
        Parameters:
            ind - int; index of filename
        returns:
            filename - string; filename of day of index
        """

        dirs = os.listdir(self.data_path)
        dirs.sort()
        filepath = self.data_path + "/" + dirs[ind]
        return(filepath)



    def multi_histmake(self, n, minfreq, maxfreq, bins_, lat1, lat0,\
                       abs = False, norm = True, pole = "north"):

        """
        make histograms of relative difference in integrated fouriers using
        multiple dates of data

        Parameters:
            n - int; number of indices in time window used in fft_time_integral_holes
            minfreq - float; lower integral limit
            maxfreq - float; higher integral limit
            bins_ - int; number of bin edges.
            lat1 - float; higher latitude limit
            lat0 - float; lower latitude limit
            abs - bool; if True, will calculate absolute difference
            norm - bool; if true, will normalize relative difference
            pole - string; "north", "south" or "both". Decides what pole we are looking at.

        returns:
            hists - list; list of histograms [BA, BC, AC]
            bins - list; center of bin values.
        """
        if norm == True:
            binlist = np.linspace(-1, 1, bins_)

        else:
            binlist = np.linspace(-5000, 5000, bins_)


        histsBA = np.zeros_like(binlist[:-1])
        histsBC = np.zeros_like(histsBA)
        histsAC = np.zeros_like(histsBC)
        for i in range(self.init_loop_index, self.end_loop_index):
            file = self.filepath(i)
            infile = MatReader(file)

            temp_hists, bins = infile.histmake(n = n, minfreq = minfreq,\
                                               maxfreq = maxfreq,\
                                               bins_ = binlist,\
                                               lat1 = lat1, lat0 = lat0,\
                                               abs = abs, norm = norm,\
                                               pole = pole)
            histsBA = histsBA + temp_hists[0]
            histsBC = histsBC + temp_hists[1]
            histsAC = histsAC + temp_hists[2]

        self.BA_shift = infile.BA_shift
        self.BC_shift = infile.BC_shift
        hists = np.array([histsBA, histsBC, histsAC])
        return(hists, bins)

    def freq_sig(self, df, jump, n, bins_, abs, norm, lat1, lat0, f1, f0, pole):
        """
        Calculates standard deviation and mean as functions of frequency
        Parameters:
            df - float; Width of frequency band
            jump - float; distance between center of frequency bands
            n - int; number of indices in time window used in fft_time_integral_holes
            f1 - float; final higher frequency limit
            f0 - float; initial lower frequency limit
            bins_ - int; number of bin edges.
            lat1 - float; higher latitude limit
            lat0 - float; lower latitude limit
            abs - bool; if True, will calculate absolute difference
            norm - bool; if true, will normalize relative difference
            pole - string; "north", "south" or "both". Decides what pole we are looking at.
        returns:
            freq0s - nparray; array containing lower integral limits
            sigmas - nparray; array on the form [sigmasBA, sigmasBC, sigmasAC]
            means - nparray; array on the form [meansBA, meansBC, meansAC]
        """
        N = int((f1 - f0)/jump)
        freq0s = np.linspace(f0, f1-df, N)
        freq1s = freq0s + df
        sigmas = np.zeros((3, N))
        means = np.zeros_like(sigmas)
        for i in range(len(freq0s)):
            print(freq0s[i]/freq0s[-1])
            minfreq = freq0s[i]
            maxfreq = freq1s[i]
            hists, bins = self.multi_histmake(n = n, minfreq = minfreq,\
                                              maxfreq = maxfreq, bins_ = bins_,\
                                              lat1 = lat1, lat0 = lat0,\
                                              abs = abs, norm = norm, pole = pole)
            for j in range(len(hists)):
                width = bins[1] - bins[0]
                hists[j] = hists[j]/np.sum(hists[j]*width)

            BA_std, BA_mean = self.pro.std_mean(hists[0], bins)
            BC_std, BC_mean = self.pro.std_mean(hists[1], bins)
            AC_std, AC_mean = self.pro.std_mean(hists[2], bins)

            sigmas[0][i] = BA_std
            sigmas[1][i] = BC_std
            sigmas[2][i] = AC_std

            means[0][i] = BA_mean
            means[1][i] = BC_mean
            means[2][i] = AC_mean


        return(freq0s, sigmas, means)










if __name__ == "__main__":
    pro = SWARMprocess()
    object = MultiMat(9, 31)
