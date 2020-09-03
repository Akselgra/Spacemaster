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


    def correlator(self, data1, data2, time, start = 0, stop = 40000,\
                   shifts = 30000, timestamps = True):
        """
        Takes 2 arrays with data and 1 array with time values.
        and produces Pearson correlation numbers for the 2 data sets,
        using different time shifts.
        Assumes data points in both data sets are gathered at the same times.

        start and stop are the indices used in slicing the data arrays.
        shifts is the number of correlation values calculated.
        """
        if timestamps == True:
            time = self.stamp_to_sec(time)#turn timestamps into seconds
        dt = time[1] - time[0] #calculates the timestep
        corr_vec = np.zeros(shifts)
        shiftvec = np.zeros(shifts)

        for i in range(len(corr_vec)):
            corr_vec[i] = pearsonr(data1[start:stop], data2[start+i:stop+i])[0]
            shiftvec[i] = dt*i #fills array with shift lengths

        return(corr_vec, shiftvec)


    def timeshift(self, data1, data2, time, start = 0, stop = 40000, shifts = 30000):
        """
        Takes 2 arrays with data and 1 array with time values.
        Calculates the most significant timeshift for the 2 data sets.
        Returns the nr of indices shifted.
        """
        corr_vec, shiftvec = self.correlator(data1, data2, time, start, stop, shifts)
        indices = np.arange(len(corr_vec))
        shift_index = indices[np.where(corr_vec == np.max(corr_vec))]
        return(shift_index[0])

    def spher_to_cart(self, lat, long, rad, deg = True):
        """
        Takes 3 arrays latitude, longitude and radius
        containing spherical coordinates and converts them
        into cartesian coordinates x, y and z
        """
        if deg == True:
            lat = lat*np.pi/180
            long = long*np.pi/180

        x = rad*np.sin(lat)*np.cos(long)
        y = rad*np.sin(lat)*np.sin(long)
        z = rad*np.cos(lat)
        return(x, y, z)

    def distance(self, x1, y1, z1, x2, y2, z2, sphere = True):
        """
        Takes coordinates of 2 bodies and finds distance between them.
        Assumes spherical coordinates if sphere is set to True.
        If spherical, x = latitude, y = longitude, z = radius
        """
        if sphere == True:
            x1, y1, z1 = self.spher_to_cart(x1, y1, z1)
            x2, y2, z2 = self.spher_to_cart(x2, y2, z2)

        x3 = x2 - x1
        y3 = y2 - y1
        z3 = z2 - z1

        r = np.sqrt(x3**2 + y3**2 + z3**2)

        return(r)
