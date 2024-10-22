from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from SWARMprocess import SWARMprocess
from time import time
from scipy.optimize import curve_fit

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

        self.shifted = False #boolean to let us assert that indices are shifted.
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

        # self.BA_shift = self.timeshift_latitude(self.latB, self.latA)
        # self.BC_shift = self.timeshift_latitude(self.latB, self.latC)
        self.BA_shift = infile["BA_shift"][0][0]
        self.BC_shift = infile["BC_shift"][0][0]

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


        self.shifted = True #changing boolean to True when function is called.

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

        assert self.shifted, "Data has not been index-shifted"

        inds = self.mlat_finder(lat1, lat0, pole)[1]
        self.inds = inds
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
    
    def comp_ind_finder(self, n, minfreq, maxfreq, lat1, lat0, abs = False, norm = True,\
                 pole = "north"):
        """
        Finds comparison indices for given parameters.
        Parameters:
            n - int; window size
            minfreq - float; lower integral limit
            maxfreq - float; higher integral limit
            abs - bool; if True, will calculate absolute difference
            norm - bool; if true, will normalize relative difference
            lat1 - float; higher latitude limit
            lat0 - float; lower altitude limit
            pole - string; "both" for both poles, "north" for north pole,
                            "south" for south pole.
        returns:
            BAdiff - array; comparison indices for BA
            BCdiff - array; comparison indices for BC
            ACdiff - array; comparison indices for AC
        """
        
        assert self.shifted, "Data has not been index-shifted"

        inds = self.mlat_finder(lat1, lat0, pole)[1]
        self.inds = inds
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
        
        return(BAdiff, BCdiff, ACdiff)

    def multi_velo_inspec(self, n = 60, lat0 = 60, lat1 = 90, pole = "north"):
        """
        Divides data set into windows of length n and finds
        along-track velocity for each window.
        Uses 50% overlap
        """
        inds = self.mlat_finder(lat1, lat0, pole)[1]
        NeA = self.NeA[inds]
        NeB = self.NeB[inds]
        NeC = self.NeC[inds]

        secondsA = self.secondsA[inds]
        secondsB = self.secondsB[inds]
        secondsC = self.secondsC[inds]


        mlatA = self.mlatA[inds]
        mlatB = self.mlatB[inds]
        mlatC = self.mlatC[inds]

        mean_range = 5
        NeA = self.meanie(NeA, mean_range)
        NeB = self.meanie(NeB, mean_range)
        NeC = self.meanie(NeC, mean_range)
        
        
        N = int((len(NeA)/n*2) - 1) #nr of windows
        
        dx = (secondsB[1]-secondsB[0])*self.velB[0]
        
        nBAs = []
        nBCs = []
        nACs = []
        
        for i in range(N):
            startind = int(i/2*n)
            stopind = int((i/2+1)*n)
            temp_NeA = NeA[startind:stopind]
            temp_NeB = NeB[startind:stopind]
            temp_NeC = NeC[startind:stopind]
            
            temp_secondsA = secondsA[startind:stopind]
            temp_secondsB = secondsB[startind:stopind]
            temp_secondsC = secondsC[startind:stopind]
            
            
            curr_timediff = np.round((temp_secondsB[1:] - temp_secondsB[:-1])-(1/self.fs))
            if np.sum(curr_timediff) > 2:
                continue
            
            gradA = (temp_NeA[1:] - temp_NeA[:-1])/dx
            gradB = (temp_NeB[1:] - temp_NeB[:-1])/dx
            gradC = (temp_NeC[1:] - temp_NeC[:-1])/dx
            
            if np.max(gradA) < 0.9:
                continue
            
            stdA = np.std(gradA)
            stdB = np.std(gradB)
            stdC = np.std(gradC)
            
            meanA = temp_secondsB[np.where(gradA == np.max(gradA))][0]
            meanB = temp_secondsB[np.where(gradB == np.max(gradB))][0]
            meanC = temp_secondsB[np.where(gradC == np.max(gradC))][0]
            
            p0A = [1, meanA, stdA]
            p0B = [1, meanB, stdB]
            p0C = [1, meanB, stdB]
            
            poptA, pcovA = curve_fit(self.gaussian, temp_secondsB[:-1], gradA, p0 = p0A)
            poptB, pcovB = curve_fit(self.gaussian, temp_secondsB[:-1], gradB, p0 = p0B)
            poptC, pcovC = curve_fit(self.gaussian, temp_secondsB[:-1], gradC, p0 = p0C)
            
            nBA = poptB[1] - poptA[1]
            nBC = poptB[1] - poptC[1]
            nAC = poptA[1] - poptC[1]
            
            nBAs.append(nBA)
            nBCs.append(nBC)
            nACs.append(nAC)
            
            
        sBA = self.BA_shift/2 #time delay BA
        sBC = self.BC_shift/2 #time delay BC
        sAC = (self.BC_shift - self.BA_shift)/2
        V = self.velA[0]
        for i in range(len(nBAs)):
            VBA = self.along_track_velo(V, sBA, nBAs[i])
            VBC = self.along_track_velo(V, sBC, nBCs[i])
            VAC = self.along_track_velo(V, sAC, nACs[i])
            
            print(VBA)
            print(VBC)
            print(VAC)
            print("________________________________________")
                
           
        
        
        
    def velo_inspec(self,ind1 = 1150*2, ind2 = 1185*2):
        """
        Inspects the velocity of a large scale structure between times
        t0 and t1.
        """
        NeA = self.NeA[ind1:ind2]
        NeB = self.NeB[ind1:ind2]
        NeC = self.NeC[ind1:ind2]

        secondsA = self.secondsA[ind1:ind2]
        secondsB = self.secondsB[ind1:ind2]
        secondsC = self.secondsC[ind1:ind2]


        mlatA = self.mlatA[ind1:ind2]
        mlatB = self.mlatB[ind1:ind2]
        mlatC = self.mlatC[ind1:ind2]

        mean_range = 5
        NeA = self.meanie(NeA, mean_range)
        NeB = self.meanie(NeB, mean_range)
        NeC = self.meanie(NeC, mean_range)

        plt.figure(0)
        plt.plot(mlatB, NeB, "r")
        plt.plot(mlatA, NeA, "g")
        plt.plot(mlatC, NeC, "b")
        plt.xlabel("Geomagnetic Latitude [Degrees]")
        plt.ylabel("Electron density [cm$^{-3}$]")
        plt.legend(["sat B", "sat A", "sat C"])
        plt.title("An interesting case")
        plt.grid("on", axis = "x")
        # plt.savefig("Figures/matfigs/interesting_case.pdf")
        plt.show()
        
        #calculate comparison index of data window
        n = len(NeA)
        NeA_ = NeA*np.hanning(n)
        NeB_ = NeB*np.hanning(n)
        NeC_ = NeC*np.hanning(n)
        fftA = np.roll(np.fft.fft(NeA)[:int(n/2)], int(n/2))
        fftB = np.roll(np.fft.fft(NeB)[:int(n/2)], int(n/2))
        fftC = np.roll(np.fft.fft(NeC)[:int(n/2)], int(n/2))
        freqs = np.linspace(-1, 1, n)[int(n/2):]
        df = freqs[1]-freqs[0]
        f = 0.1
        i = int(f/df)
        
        plt.plot(mlatB, np.hanning(n))
        plt.xlabel(" ")
        plt.ylabel("extensions")
        plt.show()
        
        plt.figure(69)
        plt.plot(np.log10(freqs[1:]), np.log10(np.abs(fftA[1:])**2), "g")
        plt.plot(np.log10(freqs[1:]), np.log10(np.abs(fftB[1:])**2), "r")
        plt.plot(np.log10(freqs[1:]), np.log10(np.abs(fftC[1:])**2), "b")
        plt.xticks([f, 0.2, 0.3])
        plt.grid("on")
        plt.xlabel("Log10 Frequency [Hz]")
        plt.ylabel("Log10 PSD")
        plt.legend(["Sat B", "Sat A", "Sat C"])
        plt.show()
        
        plt.plot(mlatB, NeB_, "r")
        plt.plot(mlatA, NeA_, "g")
        plt.plot(mlatC, NeC_, "b")
        plt.legend(["Sat B", "Sat A", "Sat C"])
        plt.show()
        
        sumA = np.sum(np.abs(fftA[i:])*df)
        sumB = np.sum(np.abs(fftB[i:])*df)
        sumC = np.sum(np.abs(fftC[i:])*df)
        
        I_BA = (sumB-sumA)/np.max([sumB, sumA])
        I_BC = (sumB-sumC)/np.max([sumB, sumC])
        I_AC = (sumA-sumC)/np.max([sumA, sumC])
        
        print("I_BA = %g" % I_BA)
        print("I_BC = %g" % I_BC)
        print("I_AC = %g" % I_AC)
        
        

        dx = (secondsB[1] - secondsB[0])*self.velB[ind1]
        der_NeA = (NeA[1:] - NeA[:-1])/dx
        der_NeB = (NeB[1:] - NeB[:-1])/dx
        der_NeC = (NeC[1:] - NeC[:-1])/dx


        #doing a gaussian fit
        from scipy.optimize import curve_fit

        def gaussian(x, amp, cen, wid):
            """1-d gaussian: gaussian(x, amp, cen, wid)"""
            return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))

        init_valsA = [4000, 1164, 0.5]
        init_valsB = [4000, 1162, 0.5]
        init_valsC = [4000, 1166, 0.5]

        poptA, pcovA = curve_fit(gaussian, secondsB[:-1],der_NeA , p0=init_valsA)
        poptB, pcovB = curve_fit(gaussian, secondsB[:-1],der_NeB , p0=init_valsB)
        poptC, pcovC = curve_fit(gaussian, secondsB[:-1],der_NeC , p0=init_valsC)


        lats = np.linspace(mlatB[0], mlatB[-1], 100)
        times = np.linspace(secondsB[0], secondsB[-1], 100)
        yA = gaussian(times, poptA[0], poptA[1], poptA[2])
        yB = gaussian(times, poptB[0], poptB[1], poptB[2])
        yC = gaussian(times, poptC[0], poptC[1], poptC[2])


        plt.figure(1)
        plt.plot(lats, yB, "--r")
        plt.plot(mlatB[:-1], der_NeB, "r")

        plt.plot(lats, yA, "--g")
        plt.plot(mlatA[:-1], der_NeA, "g")

        plt.plot(lats, yC, "--b")
        plt.plot(mlatC[:-1], der_NeC, "b")
        plt.legend(["fit B", "sat B", "fit A", "sat A", "fit C", "sat C"])

        plt.title("Density gradients")
        plt.xlabel("Geomagnetic latitude [Degrees]")
        plt.ylabel("Electron density gradient [$cm^{-3}/m$]")
        plt.grid("on", axis = "x")
        # plt.savefig("Figures/matfigs/interesting_case_deri.pdf")


        vsat = np.mean(self.velA) #velocity of satellites
        sBA = self.BA_shift/2 #time delay BA
        sBC = self.BC_shift/2 #time delay BC
        sAC = (self.BC_shift - self.BA_shift)/2

        nBA = poptB[1] - poptA[1] #time between tops in difference plot
        nBC = poptB[1] - poptC[1] #time between tops in difference plot
        nAC = poptA[1] - poptC[1]

        vBA = vsat*(1 - (sBA / (sBA - nBA))) #finding bubble velocity
        vBC = vsat*(1 - (sBC / (sBC - nBC))) #finding bubble velocity
        vAC = vsat*(1 - (sAC / (sAC - nAC)))



        NBA = int(np.abs(nBA)*2)
        NBC = int(np.abs(nBC)*2)

        fixed_NeA = self.NeA[ind1 + NBA: ind2 + NBA]
        fixed_NeB = self.NeB[ind1:ind2]
        fixed_NeC = self.NeC[ind1 + NBC: ind2 + NBC]

        fixed_NeA = self.meanie(fixed_NeA, mean_range)
        fixed_NeB = self.meanie(fixed_NeB, mean_range)
        fixed_NeC = self.meanie(fixed_NeC, mean_range)

        plt.figure(2)
        plt.plot(mlatB, fixed_NeB, "r")
        plt.plot(mlatB, fixed_NeA, "g")
        plt.plot(mlatB, fixed_NeC, "b")
        plt.xlabel("Geomagnetic Latitude [Degrees]")
        plt.ylabel("Electron density [cm$^{-3}$]")
        plt.legend(["sat B", "sat A", "sat C"])
        plt.title("An interesting case, shifted")
        plt.grid("on", axis = "x")
        # plt.savefig("Figures/matfigs/interesting_case_shifted.pdf")


        print("Bubble velocity calculated from BA = %g [m/s]" % vBA)
        print("Bubble velocity calculated from BC = %g [m/s]" % vBC)
        print("Bubble velocity calculated from AC = %g [m/s]" % vAC)
        print(nBA)
        print(nBC)
        print(nAC)
        print(sBA)
        print(sBC)
        print(sAC)





def one_period_plot():
    """
    plots one period of electron density data
    """
    file = "Data/matfiles/20131221.mat"
    object = MatReader(file)

    NeA = object.NeA
    latA = object.latA
    times = object.secondsA
    mlt = object.mltA
    ind1 = 2606 #lat inds
    ind2 = 13940 #lat inds
    
    ind1 = 3197 #mlat inds
    ind2 = 14390 #mlat inds
    
    T = ind2 - ind1
    ind1 += int(T/2)
    ind2 += int(T/2)

    latA = latA[ind1:ind2]
    NeA = NeA[ind1:ind2]
    # NeA = object.meanie(NeA, 5)
    times = times[ind1:ind2]
    mlt = mlt[ind1:ind2]
    mlt = hour_round(mlt)

    lats = np.zeros_like(latA)
    lats[0] = latA[0]
    for i in range(len(latA)-1):
        dlat = latA[i+1] - latA[i]
        if dlat < 0:
            lats[i+1] = lats[i] - dlat
        else:
            lats[i+1] = lats[i] + dlat

    lats += 90

    xticks = np.array([-90, -70, -30, 30, 70, 110, 150, 210, 250, 270]) + 90
    gridticks = np.array([-90, -70, -30, 30, 70, 77, 103, 110, 150, 210, 250, 270]) + 90
    # plt.plot(lats, NeA, ".", markersize = 1)
    # plt.plot([0, 0], [0, np.max(NeA)], "k")
    # plt.plot([30, 30], [0, np.max(NeA)], "k")
    # plt.plot([60, 60], [0, np.max(NeA)], "k")
    # plt.plot([120, 120],[0, np.max(NeA)], "k")
    # plt.plot([150, 150], [0, np.max(NeA)], "k")
    # plt.plot([167, 167], [0, np.max(NeA)], "k")
    # plt.plot([193, 193], [0, np.max(NeA)], "k")
    # plt.plot([210, 210], [0, np.max(NeA)], "k")
    # plt.plot([240, 244], [0, np.max(NeA)], "k")
    # plt.plot([300, 300], [0, np.max(NeA)], "k")
    # plt.plot([330, 330], [0, np.max(NeA)], "k")
    # plt.plot([360, 360], [0, np.max(NeA)], "k")
    # plt.xticks(xticks)
    # plt.xlabel("Geomagnetic latitude going from 0 to 360 degrees, starting and ending at south pole")
    # plt.ylabel("Electron density [cm$^{-1}$]")
    # plt.title("One SWARM satellite period")
    # plt.grid("on", axis = "x", xdata = gridticks)
    #adding letters
    x = (gridticks[:-1] + gridticks[1:])/2 - 3
    y = np.zeros_like(x) -  np.max(NeA)/40
    s = ["S", "B", "A", "B", "C", "D", "C", "B", "A", "B", "S"]
    # for i in range(len(x)):
    #     plt.text(x[i], y[i], s[i], fontsize = 10)
    # plt.savefig("Figures/swarm_period.pdf")
    # plt.show()

    # plt.plot(times, latA)
    # plt.plot(times, mlt)
    # plt.show()
    print(lats[0])
    print(lats[-1])
    
    fig, ax = plt.subplots()
    ax.plot(lats, NeA, ".", markersize = 1)
    ax.set_xticks(xticks, minor=False)
    ax.set_xticks([167, 193], minor=True)
    ax.xaxis.grid(True, which = "major")
    ax.xaxis.grid(True, which = "minor")
    for i in range(len(x)):
        ax.text(x[i], y[i], s[i], fontsize = 10)
    ax.set_xlabel("Geomagnetic latitude going from 0 to 360 degrees, starting and ending at south pole")
    ax.set_ylabel("Electron density [cm$^{-1}$]")
    ax.set_title("One Swarm satellite period")
    # plt.savefig("Figures/swarm_period.pdf")
    plt.show()
    plt.plot(mlt, NeA)
    plt.show()
    plt.plot(mlt, lats)
    plt.show()
def comparison_plotter():
    """
    Plots comparison indices
    """
    date = "20131215"
    file = "Data/matfiles/" + date + ".mat"
    object = MatReader(file)
    
    ind1 = 2606
    ind2 = int(13940 - (13940 - ind1)/2)
    
    ind1 = 0
    ind2 = -1

    n = 120
    minfreq = 0.2
    maxfreq = 1
    
    NeA = object.NeA[ind1:ind2]
    NeB = object.NeB[ind1:ind2]
    NeC = object.NeC[ind1:ind2]
    
    secondsA = object.secondsA[ind1:ind2]
    secondsB = object.secondsB[ind1:ind2]
    secondsC = object.secondsA[ind1:ind2]
    
    latA = object.latA[ind1:ind2]
    latB = object.latB[ind1:ind2]
    latC = object.latC[ind1:ind2]
    
    mlatA = object.mlatA[ind1:ind2]
    mlatB = object.mlatB[ind1:ind2]
    mlatC = object.mlatC[ind1:ind2]
    
    timesA, fftA = object.fft_time_holes_integral(NeA, secondsA, n, 2,\
                                    minfreq = minfreq, maxfreq = maxfreq)
        
    timesB, fftB = object.fft_time_holes_integral(NeB,  secondsB, n, 2,\
                                    minfreq = minfreq, maxfreq = maxfreq)
        
    timesC, fftC = object.fft_time_holes_integral(NeC, secondsC, n, 2,\
                                    minfreq = minfreq, maxfreq = maxfreq)
        
    # trans_BA = np.abs(fftB - fftA)
    # trans_BC = np.abs(fftB - fftC)
    # trans_AC = np.abs(fftA - fftC)
    
    # trans_BA = trans_BA/np.max(trans_BA)
    # trans_BC = trans_BC/np.max(trans_BC)
    # trans_AC = trans_AC/np.max(trans_AC)

    trans_BA = np.zeros_like(fftA) + 0.1
    trans_BC = np.zeros_like(fftA) + 0.1
    trans_AC = np.zeros_like(fftA) + 0.1
    
    BA_95 = np.percentile(np.abs(fftB-fftA), 95)
    BA_50 = np.percentile(np.abs(fftB-fftA), 68)
    
    BC_95 = np.percentile(np.abs(fftB-fftC), 95)
    BC_50 = np.percentile(np.abs(fftB-fftC), 68)
    
    AC_95 = np.percentile(np.abs(fftA-fftC), 95)
    AC_50 = np.percentile(np.abs(fftA-fftC), 68)
    
    print(BA_95)
    for i in range(len(trans_BA)):
        if np.abs(fftB[i] - fftA[i]) >= BA_95:
            trans_BA[i] = 1
        elif np.abs(fftB[i]-fftA[i]) >= BA_50:
            trans_BA[i] = 0.5
        
        if np.abs(fftB[i] - fftC[i]) >= BC_95:
            trans_BC[i] = 1
        elif np.abs(fftB[i]-fftC[i]) >= BC_50:
            trans_BC[i] = 0.5
            
        if np.abs(fftA[i] - fftC[i]) >= AC_95:
            trans_AC[i] = 1
        elif np.abs(fftA[i]-fftC[i]) >= AC_50:
            trans_AC[i] = 0.5
            

    # trans_BA = np.zeros_like(fftA) + 1
    # trans_BC = np.zeros_like(fftA) + 1
    # trans_AC = np.zeros_like(fftA) + 1

    comp_ind_BA = object.relative_diff(fftB, fftA, abs = False)
    comp_ind_BC = object.relative_diff(fftB, fftC, abs = False)
    comp_ind_AC = object.relative_diff(fftA, fftC, abs = False)
    
    latitudesA = object.lat_from_time(secondsA, latA, timesA)
    latitudesB = object.lat_from_time(secondsB, latB, timesB)
    latitudesC = object.lat_from_time(secondsC, latC, timesC)
    
    comp_lat_BA = (latitudesB + latitudesA)/2
    comp_lat_BC = (latitudesB + latitudesC)/2
    comp_lat_AC = (latitudesA + latitudesC)/2
    
    figs, axs = plt.subplots(3, 1, sharex = True, sharey = True)
    
    NeA = object.meanie(NeA, mean_range = 5)
    NeB = object.meanie(NeB, mean_range = 5)
    NeC = object.meanie(NeC, mean_range = 5)
    
    for i in range(len(trans_BA)):
        axs[0].plot(comp_lat_BA[i], comp_ind_BA[i], "r.", alpha = trans_BA[i])
        axs[1].plot(comp_lat_BC[i], comp_ind_BC[i], "g.", alpha = trans_BC[i])
        axs[2].plot(comp_lat_AC[i], comp_ind_AC[i], "b.", alpha = trans_AC[i])
    #axs[0].plot(comp_lat_BA, comp_ind_BA, "r.")
    #axs[1].plot(comp_lat_BC, comp_ind_BC, "g.")
    #axs[2].plot(comp_lat_AC, comp_ind_AC, "b.")
    # axs[0].plot(mlatB, NeB/np.max(NeB), "r")
    # axs[0].plot(mlatA, NeA/np.max(NeA), "g")
    # axs[1].plot(mlatB, NeB/np.max(NeB), "r")
    # axs[1].plot(mlatC, NeC/np.max(NeC), "b")
    # axs[2].plot(mlatA, NeA/np.max(NeA), "g")
    # axs[2].plot(mlatC, NeC/np.max(NeC), "b")
    
    
    
    axs[0].grid("on")
    axs[1].grid("on")
    axs[2].grid("on")
    axs[0].set_xticks([-90, -77, -70, -30, 0, 30, 70, 77, 90])
    axs[1].set_xticks([-90, -77, -70, -30, 0, 30, 70, 77, 90])
    axs[2].set_xticks([-90, -77, -70, -30, 0, 30, 70, 77, 90])
    
    xlabels = ["Latitude","Latitude","Latitude"]
    ylabels = ["$I_{BA}$","$I_{BC}$","$I_{AC}$"]
    for i in range(len(axs.flat)):
        axs.flat[i].set(xlabel=xlabels[i], ylabel=ylabels[i])
        #axs.flat[i].set_aspect("equal", "box")
        #axs.flat[i].set_xlim(-1, 1)
        axs.flat[i].set_ylim(-1, 1)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()
    
    
    
    plt.savefig("Figures/matfigs/comp_ind/comparison_indices_lat_" + date + ".pdf")
    plt.show()
    
    # timies = timesB[np.where(comp_ind_BA > 0.8)]
    # print(timies)
    # print(comp_ind_BA)
    
    # plt.plot(secondsB, NeA/np.max(NeA), "g")
    # plt.plot(secondsB, NeB/np.max(NeB), "r")
    # plt.plot(timesB, comp_ind_BA, "ko")
    # plt.show()
    

def distance_plot():
    """
    plots distance between data points
    """
    day = "20131210"
    file = "Data/matfiles/" + day + ".mat"
    object = MatReader(file)



    # ind1 = 2606
    # ind2 = 13940 + 1*7000

    ind1 = 0
    ind2 = 150000

    times = object.secondsB[ind1:ind2]

    xA = object.latA[ind1:ind2]
    yA = object.longA[ind1:ind2]
    zA = object.radA[ind1:ind2]

    xB = object.latB[ind1:ind2]
    yB = object.longB[ind1:ind2]
    zB = object.radB[ind1:ind2]

    xC = object.latC[ind1:ind2]
    yC = object.longC[ind1:ind2]
    zC = object.radC[ind1:ind2]

    mltA = object.mltA[ind1:ind2]
    mltB = object.mltB[ind1:ind2]


    dist_BA = object.great_circle_distance(xB, yB, zB, xA, yA, zA)
    dist_BC = object.great_circle_distance(xB, yB, zB, xC, yC, zC)
    dist_AC = object.great_circle_distance(xA, yA, zA, xC, yC, zC)


    plt.figure(0)
    plt.plot(times, dist_BA)
    plt.plot(times, dist_BC)
    plt.plot(times, dist_AC)
    plt.title("Distance over time")
    plt.xlabel("Seconds since midnight UTC of sat B")
    plt.ylabel("Distance [m]")
    plt.legend(["B - A", "B - C", "A - C"])
    plt.grid("on")
    # plt.savefig("Figures/matfigs/distance_over_time_" + day + ".pdf")
    plt.figure(1)

    plt.plot(times, xB - xA)
    plt.plot(times, xB - xC)
    plt.title("difference in latitude")
    plt.xlabel("time of sat B [s]")
    plt.ylabel("difference in latitude [degrees]")
    plt.legend(["B - A", "B - C"])

    plt.figure(2)

    plt.plot(times, yB - yA)
    plt.plot(times, yB - yC)
    #plt.axis([0, 7000, -10, 10])
    plt.title("difference in longitude")
    plt.xlabel("time of sat B [s]")
    plt.ylabel("difference in latitude [degrees]")
    plt.legend(["B - A", "B - C"])
    plt.show()

    # plt.plot(times, zB - zA)
    # plt.plot(times, zB - zC)
    # plt.title("Difference in altitude")
    # plt.xlabel("time of sat B [s]")
    # plt.ylabel("difference in altitude [m]")
    # plt.legend(["B - A", "B - C"])
    # plt.show()

    plt.plot(xB, dist_BA)
    plt.plot(xB, dist_BC)
    plt.plot(xB, dist_AC)
    plt.xlabel("Latitude [Degrees]")
    plt.ylabel("Distance [m]")
    plt.title("Distance over latitude")
    plt.grid("on")
    plt.legend(["B - A", "B - C", "A - C"])
    plt.xticks([-90, -77, -70, -30, 0, 30, 70, 77, 90])
    plt.savefig("Figures/matfigs/distance_over_latitude_" + day + ".pdf")
    plt.show()

    mltdiff = mltA - mltB

    for i in range(len(mltdiff)):
        if mltdiff[i] > 24:
            mltdiff[i] = mltdiff[i] - 24
        elif mltdiff[i] < -24:
            mltdiff[i] = mltdiff[i] + 24


    print("mean distance B-A = %g m" % np.mean(dist_BA))
    print("mean distance B-C = %g m" % np.mean(dist_BC))

    print("velocity times timediff BA = %g" % (object.BA_shift/2*np.mean(object.velA)))
    print("velocity times timediff BC = %g" % (object.BC_shift/2*np.mean(object.velC)))
    print(np.min(dist_BA))
    print(object.BA_shift)
    print(object.BC_shift - object.BA_shift)
    print(np.mean(object.velA))

def hour_round(hours):
    is_larger = np.nonzero(hours > 24)
    hours[is_larger] = hours[is_larger] - 24
    return(hours)

def index_shift_plot():
    """
    plots an example of index shifting
    """
    file = "Data/matfiles/20131221.mat"
    object = MatReader(file)
    thing = object.shifted == False
    assert thing, "shifter must be commented out"
    
    start = 1920 #index of 16 minutes
    stop = 3120 #index of 26 minutes
    
    NeA = object.NeA[start:stop]
    NeB = object.NeB[start:stop]
    NeC = object.NeC[start:stop]
    NeA = object.meanie(NeA, 5)
    NeB = object.meanie(NeB, 5)
    NeC = object.meanie(NeC, 5)
    
    secondsA = object.secondsA[start:stop]
    secondsB = object.secondsB[start:stop]
    secondsC = object.secondsC[start:stop]
    
    plt.plot(secondsB, NeB, "r")
    plt.plot(secondsB, NeA, "g")
    plt.plot(secondsB, NeC, "b")
    plt.xlabel("Seconds since midnight UTC of satellite B")
    plt.ylabel("Electron density [cm$^{-3}$]")
    plt.legend(["Satellite B", "Satellite A", "Satellite C"])
    plt.title("Electron density data before index-shift")
    plt.savefig("Figures/preshift_example.pdf")
    plt.show()
    
    object.shifter()
    NeA = object.NeA[start:stop]
    NeB = object.NeB[start:stop]
    NeC = object.NeC[start:stop]
    NeA = object.meanie(NeA, 5)
    NeB = object.meanie(NeB, 5)
    NeC = object.meanie(NeC, 5)
    
    secondsA = object.secondsA[start:stop]
    secondsB = object.secondsB[start:stop]
    secondsC = object.secondsC[start:stop]
    
    plt.plot(secondsB, NeB, "r")
    plt.plot(secondsB, NeA, "g")
    plt.plot(secondsB, NeC, "b")
    plt.xlabel("Seconds since midnight UTC of satellite B")
    plt.ylabel("Electron density [cm$^{-3}$]")
    plt.legend(["Satellite B", "Satellite A", "Satellite C"])
    plt.title("Electron density data after index-shift")
    plt.savefig("Figures/postshift_example.pdf")
    plt.show()
    
if __name__ == "__main__":
    fig_width_pt = 418.0  # Get this from LaTeX using \showthe\columnwidth
    # fig_width_pt = 575
    inches_per_pt = 1.0/72.27               # Convert pt to inches
    golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
    fig_width = fig_width_pt*inches_per_pt  # width in inches
    fig_height =fig_width*golden_mean       # height in inches
    fig_size = [fig_width,fig_height]
    
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Roman",
    #   "font.family": "sans-serif",
        # "font.sans-serif": ["Helvetica"],
        'font.size' : 10,
        'axes.labelsize' : 10,
        'font.size' : 10,
    #    'text.fontsize' : 10,
        'legend.fontsize': 10,
        'xtick.labelsize' : 10,
        'ytick.labelsize' : 10,
        'figure.figsize': fig_size
    })
    #matplotlib.use('pgf')

    filename = "Data/matfiles/20131221.mat"
    object = MatReader(filename)
    object.velo_inspec()
