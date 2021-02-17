from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
from time import time

class MatReader(SWARMprocess):
    """
    class to read .mat files containing SWARM data.
    """

    def __init__(self, filename):
        """
        Parameters:
            filename:
                string; name of .mat file to be read
        """

        infile = loadmat(filename)
        self.NeA = infile["NeA"][0]
        self.NeB = infile["NeB"][0]
        self.NeC = infile["NeC"][0]

        self.longA = infile["longA"][0]
        self.longB = infile["longB"][0]
        self.longC = infile["longC"][0]

        self.latA = infile["latA"][0]
        self.latB = infile["latB"][0]
        self.latC = infile["latC"][0]

        self.radA = infile["radA"][0]
        self.radB = infile["radB"][0]
        self.radC = infile["radC"][0]

        self.velA = infile["velA"][0]
        self.velB = infile["velB"][0]
        self.velC = infile["velC"][0]

        self.altA = infile["altA"][0]
        self.altB = infile["altB"][0]
        self.altC = infile["altC"][0]

        self.mlatA = infile["mlatA"][0]
        self.mlatB = infile["mlatB"][0]
        self.mlatC = infile["mlatC"][0]

        self.mlongA = infile["mlongA"][0]
        self.mlongB = infile["mlongB"][0]
        self.mlongC = infile["mlongC"][0]

        self.mltA = infile["mltA"][0]
        self.mltB = infile["mltB"][0]
        self.mltC = infile["mltC"][0]

        self.secondsA = infile["secondsA"][0]
        self.secondsB = infile["secondsB"][0]
        self.secondsC = infile["secondsC"][0]

        self.fs = 2

        self.BA_shift = self.timeshift_latitude(self.latB, self.latA)
        self.BC_shift = self.timeshift_latitude(self.latB, self.latC)

        self.shifter()

    def shifter(self):
        """
        shifts data so as to minimize distance in latitude.
        Lets B be stationary while index shifting A and C forward
        so that
        A[t] -> A[t - t_BA]
        B[t] -> B[t]
        C[t] -> C[t - t_BC]
        """
        #self.BA_shift = self.timeshift_latitude(self.latB, self.latA)
        #self.BC_shift = self.timeshift_latitude(self.latB, self.latC)




        secondsA = self.secondsA
        secondsB = self.secondsB
        secondsC = self.secondsC

        NeA = self.holefill(self.NeA, secondsA)
        NeB = self.holefill(self.NeB, secondsB)
        NeC = self.holefill(self.NeC, secondsC)

        start = 0
        stop = len(NeA) - np.max(np.array([self.BA_shift, self.BC_shift]))

        startA = start + self.BA_shift
        stopA = stop + self.BA_shift

        startC = start + self.BC_shift
        stopC = stop + self.BC_shift

        NeA = NeA[startA:stopA]
        NeB = NeB[start:stop]
        NeC = NeC[startC:stopC]

        longA = self.holefill(self.longA, secondsA)
        longB = self.holefill(self.longB, secondsB)
        longC = self.holefill(self.longC, secondsC)
        longA = longA[startA:stopA]
        longB = longB[start:stop]
        longC = longC[startC:stopC]

        latA = self.holefill(self.latA, secondsA)
        latB = self.holefill(self.latB, secondsB)
        latC = self.holefill(self.latC, secondsC)
        latA = latA[startA:stopA]
        latB = latB[start:stop]
        latC = latC[startC:stopC]

        radA = self.holefill(self.radA, secondsA)
        radB = self.holefill(self.radB, secondsB)
        radC = self.holefill(self.radC, secondsC)
        radA = radA[startA:stopA]
        radB = radB[start:stop]
        radC = radC[startC:stopC]

        velA = self.holefill(self.velA, secondsA)
        velB = self.holefill(self.velB, secondsB)
        velC = self.holefill(self.velC, secondsC)
        velA = velA[startA:stopA]
        velB = velB[start:stop]
        velC = velC[start:stop]

        altA = self.holefill(self.altA, secondsA)
        altB = self.holefill(self.altB, secondsB)
        altC = self.holefill(self.altC, secondsC)
        altA = altA[startA:stopA]
        altB = altB[start:stop]
        altC = altC[startC:stopC]


        mlatA = self.holefill(self.mlatA, secondsA)
        mlatB = self.holefill(self.mlatB, secondsB)
        mlatC = self.holefill(self.mlatC, secondsC)
        mlatA = mlatA[startA:stopA]
        mlatB = mlatB[start:stop]
        mlatC = mlatC[startC:stopC]

        mlongA = self.holefill(self.mlongA, secondsA)
        mlongB = self.holefill(self.mlongB, secondsB)
        mlongC = self.holefill(self.mlongC, secondsC)
        mlongA = mlongA[startA:stopA]
        mlongB = mlongB[start:stop]
        mlongC = mlongC[startC:stopC]

        mltA = self.holefill(self.mltA, secondsA)
        mltB = self.holefill(self.mltB, secondsB)
        mltC = self.holefill(self.mltC, secondsC)
        mltA = mltA[startA:stopA]
        mltB = mltB[start:stop]
        mltC = mltC[startC:stopC]

        secondsA = self.holefill(secondsA, secondsA)
        secondsB = self.holefill(secondsB, secondsB)
        secondsC = self.holefill(secondsC, secondsC)
        secondsA = secondsA[startA:stopA]
        secondsB = secondsB[start:stop]
        secondsC = secondsC[startC:stopC]

        indsA = np.nonzero(secondsA)[0]
        indsB = np.nonzero(secondsB)[0]
        indsC = np.nonzero(secondsC)[0]

        inds = np.intersect1d(indsA, indsB)
        inds = np.intersect1d(inds, indsC)

        self.NeA = NeA[inds]
        self.NeB = NeB[inds]
        self.NeC = NeC[inds]

        self.longA = longA[inds]
        self.longB = longB[inds]
        self.longC = longC[inds]

        self.latA = latA[inds]
        self.latB = latB[inds]
        self.latC = latC[inds]

        self.radA = radA[inds]
        self.radB = radB[inds]
        self.radC = radC[inds]

        self.velA = velA[inds]
        self.velB = velB[inds]
        self.velC = velC[inds]

        self.altA = altA[inds]
        self.altB = altB[inds]
        self.altC = altC[inds]

        self.mlatA = mlatA[inds]
        self.mlatB = mlatB[inds]
        self.mlatC = mlatC[inds]

        self.mlongA = mlongA[inds]
        self.mlongB = mlongB[inds]
        self.mlongC = mlongC[inds]

        self.mltA = mltA[inds]
        self.mltB = mltB[inds]
        self.mltC = mltC[inds]

        self.secondsA = secondsA[inds]
        self.secondsB = secondsB[inds]
        self.secondsC = secondsC[inds]

    def mlat_finder(self, lat1, lat0, pole = "north"):
            """
            Splits data set in high and low latitude regions
            split by geomagnetic latitude
            returns:
                indsA - list; list on the form [A0, A1]
                               where A0...C1 are arrays of indices
                               with 1 being the region between lat0 and lat1
                               and 0 being the region between lat1 and lat0
            """

            if pole == "both":
                lowerA = np.abs(self.mlatA) < lat1
                higherA = np.abs(self.mlatA) > lat0
                is_poleA = lowerA * higherA

            elif pole == "north":
                lowerA = (self.mlatA) < lat1
                higherA = (self.mlatA) > lat0
                is_poleA = lowerA * higherA

            elif pole == "south":
                lowerA = (self.mlatA) > lat1
                higherA = (self.mlatA) < lat0
                is_poleA = lowerA * higherA

            high_lat_A = np.where(is_poleA == 1)
            low_lat_A = np.where(is_poleA == 0)
            indsA = [low_lat_A, high_lat_A]

            return indsA


    def histmake(self, n, minfreq, maxfreq, bins_, lat1, lat0, abs = False, norm = True,\
                 pole = "north"):
        """
        Makes histogram of integrated fourier spectrum for the given
        latitude region.
        Parameters:
            n - int; window size
            minfreq - float; lower integral limit
            maxfreq - float; higher integral limit
            bins_ - array or int; bins to be given to numpy.histogram
            abs - bool; if True, will calculate absolute difference
            norm - bool; if true, will normalize relative difference
            lat1 - float; higher latitude limit
            lat0 - float; lower altitude limit
            pole - string; "both" for both poles, "north" for north pole,
                            "south" for south pole.

        returns:
            hists - list; list of histograms on the form [BA, BC, AC]
            bins - list; center of bin values.
        """

        inds = self.mlat_finder(lat1, lat0, pole)[1]
        NeA = self.NeA[inds]
        NeB = self.NeB[inds]
        NeC = self.NeC[inds]
        secondsA = self.secondsA[inds]
        secondsB = self.secondsB[inds]
        secondsC = self.secondsC[inds]

        timesA, fftA = self.fft_time_holes_integral(signal = NeA,\
                       seconds = secondsA, n = n, fs = self.fs, minfreq\
                       = minfreq, maxfreq = maxfreq)

        timesB, fftB = self.fft_time_holes_integral(signal = NeB,\
                       seconds = secondsB, n = n, fs = self.fs, minfreq\
                       = minfreq, maxfreq = maxfreq)

        timesC, fftC = self.fft_time_holes_integral(signal = NeC,\
                       seconds = secondsC, n = n, fs = self.fs, minfreq\
                       = minfreq, maxfreq = maxfreq)

        BAdiff = self.relative_diff(fftB, fftA, norm = norm, abs = abs)
        BCdiff = self.relative_diff(fftB, fftC, norm = norm, abs = abs)
        ACdiff = self.relative_diff(fftA, fftC, norm = norm, abs = abs)

        BAhist, bins = np.histogram(BAdiff, bins = bins_)
        BChist, bins = np.histogram(BCdiff, bins = bins_)
        AChist, bins = np.histogram(ACdiff, bins = bins_)

        bins = (bins[1:] + bins[:-1])/2
        hists = [BAhist, BChist, AChist]

        return(hists, bins)


if __name__ == "__main__":

    def hour_round(hours):
        is_larger = np.nonzero(hours > 24)
        hours[is_larger] = hours[is_larger] - 24
        return(hours)
    file = "Data/matfiles/20131214.mat"
    object = MatReader(file)
    inds = object.mlat_finder(75, 65)
    NeA = object.NeA[inds[1]]
    secondsA = object.secondsA[inds[1]]
    times, fourier_int = object.fft_time_holes_integral(NeA, secondsA, 100, 2, minfreq = 0.1, maxfreq = 0.9)

    plt.figure(0)
    plt.plot(times/60/60, fourier_int, ".")
    plt.xlabel("UTC")
    plt.ylabel("integrated fourier")

    plt.figure(1)
    plt.xlabel("MLT")
    plt.ylabel("Ne")
    plt.plot(object.mltA[inds[1]], NeA, ".")

    def hour(x):
        inds = np.nonzero(x > 24)
        x[inds] -= 24
        return x

    plt.figure(2)
    plt.plot(object.mlatA, hour(object.mltA + object.secondsA/60/60))
    plt.xlabel("Geomagnetic Latitude")
    plt.ylabel("Magnetic Local Time")
    plt.show()
    """
    hists, bins = object.histmake(100, 0.1, 0.9, np.linspace(-1, 1, 50),75,65)
    BAhist = hists[0]
    width = bins[1] - bins[0]
    plt.bar(bins, BAhist, width = width)
    plt.show()
    """
