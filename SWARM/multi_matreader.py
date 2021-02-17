from SWARMprocess import SWARMprocess
import numpy as np
import matplotlib.pyplot as plt
from general_matreader import MatReader
import os
from scipy.io import loadmat



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
        
        self.data_path = data_path #path to .mat files
        
        self.init_loop_index = day0 - day_0
        self.end_loop_index = day1 - day_0 + 1
        
        
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
            n - int; number of indices in time window used in fft_time_integral
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
            binlist = np.linspace(-250000, 250000, bins_)
        
        
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
        
        hists = np.array([histsBA, histsBC, histsAC])
        return(hists, bins)
        
        
        
        



if __name__ == "__main__":
    pro = SWARMprocess()
    object = MultiMat(15, 31)
    hists, bins = object.multi_histmake(100, 0.2, 0.4, 50, 75, 65)
    width = bins[1] - bins[0]
    hists[0] = hists[0]/np.sum(hists[0]*width)
    std, mean = pro.std_mean(hists[0], bins)
    xs = np.linspace(-1, 1, 1000)
    plt.bar(bins, hists[0], width = width)
    plt.plot(xs, pro.gauss_curve(xs, mean, std), "r")
    plt.show()