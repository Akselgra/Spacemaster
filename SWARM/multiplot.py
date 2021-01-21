import multi_SWARMreader
import numpy as np
import matplotlib.pyplot as plt
import time

def sigma_plotter():
    start = time.time()
    object = multi_SWARMreader.MultiSWARM(2013, 12, 9, 2013, 12, 31)
    freq0s, sigmas, means = object.freq_sig(df = 0.01, n = 100, bins_ = 50)

    plt.plot(freq0s, sigmas[0])
    plt.plot(freq0s, sigmas[1])
    plt.plot(freq0s, sigmas[2])
    plt.xlabel("Lower integral limit")
    plt.ylabel("STD")
    plt.legend(["B-A, B-C, A-C"])
    plt.title("standard deviations of histograms, df = 0.01")
    plt.show()
    stop = time.time()
    print(stop-start)


def sigma_plotter_lat():
    df = 0.05
    start = time.time()
    object = multi_SWARMreader.MultiSWARM(2013, 12, 25, 2013, 12, 25)
    freq0s, sigmas, means = object.freq_sig_lat(df = df, n = 100, bins_ = 50, lat1 = 90, lat0 = 75)

    plt.plot(freq0s, sigmas[0])
    plt.plot(freq0s, sigmas[1])
    plt.plot(freq0s, sigmas[2])
    plt.plot(freq0s, sigmas[3])
    plt.plot(freq0s, sigmas[4])
    plt.plot(freq0s, sigmas[5])
    plt.xlabel("Lower integral limit")
    plt.ylabel("STD")
    plt.legend(["BAH", "BAL", "BCH", "BCL", "ACH", "ACL"])
    plt.title("standard deviations of histograms, df = %g, 25.12.13" % df)
    plt.show()
    stop = time.time()
    print(stop-start)

def spec_lat_hist_plotter():
    start = time.time()
    object = multi_SWARMreader.MultiSWARM(2013, 12, 24, 2013, 12, 24)
    minfreq = 0
    maxfreq = 1
    hists, bins = object.multi_histmake_lat(bins_ = 50, minfreq = minfreq,\
                                            maxfreq = maxfreq, lat0 = 75, lat1 = 90,\
                                            norm = True, n = 100)



    means = np.zeros(len(bins))
    stds = np.zeros_like(means)
    for i in range(len(bins)):
        dbins = bins[i][1]-bins[i][0]
        hists[i] = hists[i]/np.sum(hists[i]*dbins)
        curr_std, curr_mean = object.pro.std_mean(hists[i], bins[i])
        means[i] = curr_mean
        stds[i] = curr_std


    xs = np.linspace(np.min(bins[0]), np.max(bins[0]), 1000)
    for i in range(len(hists)):
        gauss = object.pro.gauss_curve(x = xs, mean = means[i], std = stds[i])
        plt.figure(i)
        width = bins[i][1] - bins[i][0]
        plt.plot(xs, gauss, "r")
        plt.bar(bins[i], hists[i], width = 0.9*width)
        plt.xlabel("Relative difference")
        plt.ylabel("Normalized occurence")
        plt.legend(["Gauss approximation", "Histogram"])
        stringies = ["B-A high", "B-A low", "B-C high", "B-C low", "A-C high", "A-C low"]
        plt.title(stringies[i])

    plt.show()

    BA_shift = object.BA_shift
    BC_shift = object.BC_shift
    AC_shift = BC_shift - BA_shift
    BA_shift = BA_shift/2
    BC_shift = BC_shift/2
    AC_shift = AC_shift/2
    plt.plot([BA_shift, AC_shift, BC_shift], [stds[0], stds[4], stds[2]], "o-")
    plt.xlabel("Time difference between satellites")
    plt.ylabel("standard deviation in high latitude histograms")
    plt.show()


    stop = time.time()
    print(stop-start)

sigma_plotter_lat()
