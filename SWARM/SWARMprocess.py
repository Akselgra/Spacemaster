import numpy as np
from scipy.stats import pearsonr


class SWARMprocess():
    """
    Class for processing SWARM data
    """

    def __init__(self):
        self.x = 5


    def stamp_to_sec(self, times):
        """
        Takes timestamps and returns seconds.
        """
        seconds = np.zeros(len(times))

        for i in range(len(times)):
            curr_str = str(times[i])
            hrs = float(curr_str[11:13])
            mins = float(curr_str[14:16])
            secs = float(curr_str[17:])
            seconds[i] = 3600*hrs + 60*mins + secs

        return(seconds)


    def correlator(self, data1, data2, time, start = 0, stop = 40000, shifts = 30000):
        """
        Takes 2 arrays with data and 1 array with time values.
        and produces Pearson correlation numbers for the 2 data sets,
        using different time shifts.
        Assumes data points in both data sets are gathered at the same times.

        start and stop are the indices used in slicing the data arrays.
        shifts is the number of correlation values calculated.
        """

        time = self.stamp_to_sec(time)#turn timestamps into seconds
        dt = time[1] - time[0] #calculates the timestep
        corr_vec = np.zeros(shifts)
        shiftvec = np.zeros(shifts)

        for i in range(len(corr_vec)):
            corr_vec[i] = pearsonr(data1[start:stop], data2[start+i:stop+i])[0]
            shiftvec[i] = dt*i #fills array with shift lengths

        return(corr_vec, shiftvec)
