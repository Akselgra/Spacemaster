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
    df = 0.01
    start = time.time()
    object = multi_SWARMreader.MultiSWARM(2013, 12, 9, 2013, 12, 31)
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
    plt.title("standard deviations of histograms, df = %g" % df)
    plt.show()
    stop = time.time()
    print(stop-start)


sigma_plotter_lat()
