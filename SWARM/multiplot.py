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
    plt.show()
    stop = time.time()
    print(stop-start)


sigma_plotter()
