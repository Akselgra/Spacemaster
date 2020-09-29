import numpy as np
from scipy.stats import pearsonr


class SWARMprocess():
    """
    Class for processing SWARM data
    """

    def __init__(self):
        self.Re = 6371008.8 #earth radius


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


    def correlator(self, data1, data2, time = False, start = 0, stop = 40000,\
                   shifts = 30000, timestamps = True):
        """
        Takes 2 arrays with data and 1 array with time values.
        and produces Pearson correlation numbers for the 2 data sets,
        using different time shifts.
        Assumes data points in both data sets are gathered at the same times.

        start and stop are the indices used in slicing the data arrays.
        shifts is the number of correlation values calculated.

        If time == False: returns indices shifted instead of seconds.
        """
        if time != False:
            if timestamps == True:
                time = self.stamp_to_sec(time)#turn timestamps into seconds
            dt = time[1] - time[0] #calculates the timestep
            corr_vec = np.zeros(shifts)
            shiftvec = np.zeros(shifts)

            for i in range(len(corr_vec)):
                corr_vec[i] = pearsonr(data1[start:stop], data2[start+i:stop+i])[0]
                shiftvec[i] = dt*i #fills array with shift lengths

            return(corr_vec, shiftvec)

        else:
            corr_vec = np.zeros(shifts)
            shiftvec = np.zeros(shifts)

            for i in range(len(corr_vec)):
                corr_vec[i] = pearsonr(data1[start:stop], data2[start+i:stop+i])[0]
                shiftvec[i] = i
            return(corr_vec, shiftvec)

    def timeshift(self, data1, data2, time = False, start = 0, stop = 40000, shifts = 30000):
        """
        Takes 2 arrays with data and 1 array with time values.
        Calculates the most significant timeshift for the 2 data sets.
        Returns the nr of indices shifted.
        """
        corr_vec, shiftvec = self.correlator(data1, data2, time, start, stop, shifts)
        indices = np.arange(len(corr_vec))
        shift_index = indices[np.where(corr_vec == np.max(corr_vec))]
        return(shift_index[0])

    def timeshift_latitude(self, lat1, lat2, start = 0, stop = 40000, shifts = 30000):
        """
        Assuming only latitude changes, takes 2 arrays with latitudes.
        Calculates the index shift that gives the smallest difference in latitudes
        With a change in only latitudes, this gives the measurement points with
        the smallest change in distance.
        """
        meandist = np.zeros(shifts)
        indices = np.arange(shifts)

        for i in range(shifts):
            meandist[i] = np.mean(np.abs(lat1[start:stop] - lat2[start + i: stop + i]))

        shift_index = indices[np.where(meandist == np.min(meandist))]
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

    def distance(self, x1, y1, z1, x2, y2, z2, sphere = True, rad = False):
        """
        Takes coordinates of 2 bodies and finds distance between them.
        Assumes spherical coordinates if sphere is set to True.
        If spherical, x = latitude, y = longitude, z = radius
        Assumes x1, ..., z2 are arrays
        """
        if sphere == True: #distance in spherical coordinates
            if rad == False:
                x1 = x1*np.pi/180
                y1 = y1*np.pi/180
                x2 = x2*np.pi/180
                y2 = y2*np.pi/180

            a = z1**2 + z2**2
            b = np.sin(x1)*np.sin(x2)*np.cos(y1)*np.cos(y2)
            c = np.sin(x1)*np.sin(x2)*np.sin(y1)*np.sin(y2) + np.cos(x1)*np.cos(x2)
            return(np.sqrt(a - 2*z1*z2*(b + c)))

        x3 = x2 - x1
        y3 = y2 - y1
        z3 = z2 - z1

        r = np.sqrt(x3**2 + y3**2 + z3**2)

        return(r)

    def reshape(self, array, shape):
        """
        Takes 1d array and shape
        reshapes array to be of given shape.
        Array given must be 1d.
        """

        new_array = np.zeros(shape)
        if len(array) < shape[0]*shape[1]:
            app_array = np.zeros((shape[0]*shape[1] - len(array)))
            array = np.concatenate((array, app_array))

        for i in range(len(new_array)):
            for j in range(len(new_array[0])):
                new_array[i][j] = array[shape[1]*i + j]

        return(new_array)

    def meanie(self, array, mean_range):
        """
        Takes an array and a mean_range.
        Returns an array where every element is the average of the 2*mean_range
        elements neighbouring it and itself.
        """
        new_array = np.zeros(np.shape(array))
        #setting the boundaries
        for i in range(mean_range):
            new_array[i] = array[i]
        for i in range(1, mean_range+1):
            new_array[-i] = array[-i]

        for i in range(mean_range, len(array)-(mean_range)):
            new_array[i] = np.mean(array[i - mean_range:i + mean_range+1])

        return(new_array)

if __name__ == "__main__":
    pass
