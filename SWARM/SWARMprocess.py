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
        if type(time) != type(False):
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

    def correlator2(self, data1, data2):
        """
        Takes data1 and data2 and calculates pearson R correlation coefficients
        for different timeshifts.

        Assumes that data1 and data2 are of equal sizes.
        Uses data1 as background while shifting data2.

        Assumes that data1[0] and data2[0] are measured at the same times.

        Returns corrvec; array of same size as data which contains
                         pearson r correlation coefficients for
                         different shifts, going from
                         -len(data1)/2 to len(data1)/2
        """

        if len(data1) != len(data2):
            raise ValueError("data1 and data2 are not of equal length")

        n = int(len(data1)/2)
        corrvec = np.zeros(2*n)


        for i in range(n):
            corrvec[i] = pearsonr(data1[i:n+i], data2[n:2*n])[0]
            corrvec[n+i] = pearsonr(data1[i:n+i], data2[:n])[0]

        return corrvec


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

    def timeshift_latitude(self, lat1, lat2, start = 0, stop = 40000, shifts = 7000):
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
            new_array[i] = 0
        for i in range(1, mean_range+1):
            new_array[-i] = 0

        for i in range(mean_range, len(array)-(mean_range)):
            new_array[i] = np.mean(array[i - mean_range:i + mean_range+1])

        return(new_array)


    def maxdiff(self, array1, array2):
        """
        Takes 2 1d arrays of same length.
        Returns number of indices between max values in arrays
        on the form (max1, (max1 - max2)), as well as index of max value of array1.
        """
        #assert arrays equal length
        testy = len(array1) == len(array2)
        msg = "Arrays are not of equal length"
        assert testy, msg

        indices = np.arange(len(array1))
        max1 = indices[np.where(array1 == np.max(array1))]
        max2 = indices[np.where(array2 == np.max(array2))]

        max1 = np.mean(max1)
        max2 = np.mean(max2)

        return(int(max1), int(max1 - max2))

    def wavefront_finder(self, array1, array2, array3, mean_range = 10,\
        partsize = 50):
        """
        Takes 3 syncronized data sets of same length and
        mean_range, the number of data points used to smooth out the gradient
        partsize, the number of data points inspected at any one time.
        Returns:
        index of wavefronts in array 1,
        distance (indices) between wavefronts in array 1 and 2
        distance (indices) between wavefronts in array 1 and 3
        """

        #calculate gradients
        diff1 = array1[1:] - array1[:-1]
        diff2 = array2[1:] - array2[:-1]
        diff3 = array3[1:] - array3[:-1]
        #smoothen the gradients
        mean1 = self.meanie(diff1, mean_range)
        mean2 = self.meanie(diff2, mean_range)
        mean3 = self.meanie(diff3, mean_range)

        #cutting the arrays into parts of length partsize, with 50% overlap
        m = len(mean1)
        parts = int(m/partsize*2) -1

        diff_inds = []
        diffs_1_2 = []
        diffs_1_3 = []

        for i in range(parts):
            startind = int(i/2*partsize)
            stopind = int((i/2+1)*partsize)
            curr1 = mean1[startind:stopind]
            curr2 = mean2[startind:stopind]
            curr3 = mean3[startind:stopind]
            #finding distance between max values along with index of max values
            ind_1, diff_1_2 = self.maxdiff(curr1, curr2)
            ind_1 += startind
            ind_2, diff_1_3 = self.maxdiff(curr1, curr3)
            diff_inds.append(ind_1)
            diffs_1_2.append(diff_1_2)
            diffs_1_3.append(diff_1_3)

        diff_inds = np.array(diff_inds) +1
        diffs_1_2 = np.array(diffs_1_2)
        diffs_1_3 = np.array(diffs_1_3)

        return(diff_inds, diffs_1_2, diffs_1_3)



    def gauss_curve(self, x, mean, std):
        """
        Gaussian bell function
        """
        first = 1/(std*np.sqrt(2*np.pi))
        second = -0.5*(x - mean)**2/(std**2)
        return(first*np.exp(second))

    def cross_spectrum(self, u, v, fs = 2):
        """
        Takes 2 signals u and v
        returns cross spectrum of u and v
        DTFT[u]* times DTFT[v]
        """
        N = len(u)
        fftu = np.fft.fft(u)/N
        fftv = np.fft.fft(v)/N

        return(2*N/fs*np.conjugate(fftu)*fftv)


    def time_displace(self, phase, freq):
        """
        Takes phase and frequency from cross spectrum
        returns time displacement
        """

        return(phase/(2*np.pi*freq))

    def fft_time(self, signal, n, fs):
        """
        Splits signal into n windows and calculates fourier transform of each
        uses 50% overlap
        parameters:
            signal - signal to be processed
            n - width of window
            fs - sampling frequency
        returns:
            Freqs - 2D array with frequency coordinates
            Times - 2D array with times coordinates
            ffts - array containing fourier transforms
        """

        N = len(signal)

        m = 0
        n_temp = 0
        while n_temp + int(n/2) <= N:
            m += 1
            n_temp += int(n/2)

        ffts = []
        for i in range(1, m):
            ind1 = int(n/2)*(i - 1)
            ind2 = int(n/2)*(i + 1)
            curr_dat = signal[ind1:ind2]
            curr_dat = curr_dat*np.hanning(len(curr_dat))
            ffts.append(np.fft.fft(curr_dat)[:int(n/2)])

        ffts = np.array(ffts)

        times = np.arange(m-1)*int(n/2)/fs + int(n/2)/fs
        freqs = np.linspace(-fs/2, fs/2, n)[int(n/2):]
        Freqs, Times = np.meshgrid(freqs, times)

        return(Freqs, Times, ffts)

    def CSD_time(self, u, v, n, fs):
        """
        Splits signals into n windows and calculates cross spectrum of each
        uses 50% overlap
        parameters:
            u - signal1 to be processed
            v - signal2 to be processed
            n - width of window
            fs - sampling frequency
        returns:
            Freqs - 2D array with frequency coordinates
            Times - 2D array with times coordinates
            CSDs - array containing CSDs
        """

        N = len(u)

        m = 0
        n_temp = 0
        while n_temp + int(n/2) <= N:
            m += 1
            n_temp += int(n/2)

        CSDs = []
        for i in range(1, m):
            ind1 = int(n/2)*(i - 1)
            ind2 = int(n/2)*(i + 1)
            curr_dat_u = u[ind1:ind2]
            curr_dat_v = v[ind1:ind2]
            curr_dat_u = curr_dat_u*np.hanning(len(curr_dat_u))
            curr_dat_v = curr_dat_v*np.hanning(len(curr_dat_v))
            CSDs.append(self.cross_spectrum(curr_dat_u, curr_dat_v, fs)[:int(n/2)])

        CSDs = np.array(CSDs)

        times = np.arange(m-1)*int(n/2)/fs + int(n/2)/fs
        freqs = np.linspace(-fs/2, fs/2, n)[int(n/2):]
        Freqs, Times = np.meshgrid(freqs, times)

        return(Freqs, Times, CSDs)

    def fft_time_integral(self, signal, n, fs, minfreq = 0, maxfreq = 1):
        """
        Calls fft_time and integrates results over frequency.
        parameters:
            signal - signal to be processed
            n - width of window
            fs - sampling frequency
            minfreq - min integral limit
            maxfreq - max integral limit
        returns:
            times - array of time points
            fourier_int - array with integrated fourier values
        """

        assert maxfreq<=(fs/2), "maxfreq must be lower or equal nyquist frequency"

        N = len(signal)

        Freqs, Times, ffts = self.fft_time(signal, n, fs)

        ffts = np.abs(ffts) + 1e-16
        freqs = Freqs[0, :]
        df = freqs[1] - freqs[0]
        times = Times[:, 0]

        N_maxfreq = int(maxfreq/df)
        N_minfreq = int(minfreq/df)

        temp_ffts = ffts[:int(len(times)), N_minfreq:N_maxfreq]
        fourier_int = np.sum(temp_ffts, axis = 1)*df

        return(times, fourier_int)

    def CSD_time_integral(self, u, v, n, fs, minfreq = 0, maxfreq = 1):
        """
        Calls fft_time and integrates results over frequency.
        parameters:
            u - signal1 to be processed
            v - signal2 to be processed
            n - width of window
            fs - sampling frequency
            minfreq - min integral limit
            maxfreq - max integral limit
        returns:
            times - array of time points
            CSD_int - array with integrated CSDs
        """

        assert maxfreq<=(fs/2), "maxfreq must be lower or equal nyquist frequency"

        N = len(u)

        Freqs, Times, CSDs = self.CSD_time(u, v, n, fs)

        CSDs = np.abs(CSDs) + 1e-16
        freqs = Freqs[0, :]
        df = freqs[1] - freqs[0]
        times = Times[:, 0]

        N_maxfreq = int(maxfreq/df)
        N_minfreq = int(minfreq/df)

        temp_CSDs = CSDs[:int(len(times)), N_minfreq:N_maxfreq]
        CSD_int = np.sum(temp_CSDs, axis = 1)*df

        return(times, CSD_int)

    def relative_diff(self, array1, array2):
        """
        Takes 2 arrays and returns relative difference
        |array1 - array2|/max(array1, array2)

        parameters:
            array1 - array
            array2 - array
        returns:
            diff - array
        """

        assert len(array1)==len(array2), "array1 and array2 must be of equal length"

        diff = np.abs(array1 - array2)

        normals = np.zeros_like(array1) #normalization array

        for i in range(len(normals)):
            normals[i] = np.max(np.array([array1[i], array2[i]]))

        diff = diff/normals
        return(diff)


if __name__ == "__main__":
    pass
