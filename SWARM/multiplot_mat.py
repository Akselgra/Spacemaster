import multi_matreader
import SWARMprocess
import numpy as np
import matplotlib.pyplot as plt
import time
pro = SWARMprocess.SWARMprocess()


def std_timeshift_mat():
    start = time.time()
    minfreq = 0.1
    maxfreq = 0.3
    day0 = 9
    day1 = 31
    lat1 = 90
    lat0 = 77
    n = 100
    shift_list = []
    std_list = []
    k_high = []
    
    for day in range(day0, day1+1):
        object = multi_matreader.MultiMat(day, day)
        hists, bins = object.multi_histmake(n = n, minfreq = minfreq, maxfreq = maxfreq, bins_ = 50, lat1 = lat1, lat0 = lat0)
        means = np.zeros(len(hists))
        stds = np.zeros_like(means)
        for i in range(len(hists)):
            width = bins[1] - bins[0]
            hists[i] = hists[i]/np.sum(hists[i]*width)
            curr_std, curr_mean = pro.std_mean(hists[i], bins)
            means[i] = curr_mean
            stds[i] = curr_std
            
        BA_shift = object.BA_shift
        BC_shift = object.BC_shift
        AC_shift = BC_shift - BA_shift
        BA_shift = BA_shift/2
        BC_shift = BC_shift/2
        AC_shift = AC_shift/2
        shift_list.append(BA_shift)
        std_list.append(stds[0])
        shift_list.append(BC_shift)
        std_list.append(stds[1])
        shift_list.append(AC_shift)
        std_list.append(stds[2])
        if day == 14 or day == 25:
            for j in range(3):
                k_high.append(1)
        else:
            for j in range(3):
                k_high.append(0)
        
        
    
    p = np.polyfit(shift_list, std_list, deg = 1)
    a = p[0]; b = p[1]
    print("Slope of regression is %g" % a)
    print("Constant of linear regression is %g" % b)
    
    print(a/b*100)
    xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
    plt.figure(1)
    plt.plot(xs, xs*a + b)
    plt.plot(shift_list, std_list, "ko")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("Standard deviation of histograms")
    plt.title("Stds per day, f: %g - %g, a = %g, b = %g" % (minfreq, maxfreq, a, b))
    plt.legend(["Linear regression", "data points"])
    
    plt.figure(2)
    for i in range(len(shift_list)):
        if k_high[i] == 0:
            plt.plot(shift_list[i], std_list[i], "ko")
        else:
            plt.plot(shift_list[i], std_list[i], "ro")


    plt.xlabel("Time between satellites [s]")
    plt.ylabel("Standard deviation of histograms")
    plt.title("Stds per day, integral limits: %g to %g" % (minfreq, maxfreq))
    plt.show()
    
    
    
    stop = time.time()
    print(stop-start)
    

def sigma_plotter_mat():
    """
    plots std as a function of frequency
    """
    f0 = 0
    f1 = 1
    n = 100
    df = 0.1
    jump = 0.1
    lat1 = 77
    lat0 = 65
    start = time.time()
    object = multi_matreader.MultiMat(9, 31)
    freq0s, sigmas, means = object.freq_sig(df = df, jump = jump, n = n,\
                                            f0 = f0, f1 = f1, bins_ = 50,\
                                            abs = False, norm = True,\
                                            lat1 = lat1, lat0 = lat0,\
                                                pole = "north")

    plt.plot(freq0s, sigmas[0])
    plt.plot(freq0s, sigmas[1])
    plt.plot(freq0s, sigmas[2])
    plt.xlabel("Lower integral limit")
    plt.ylabel("STD")
    plt.legend(["B-A", "B-C", "A-C"])
    plt.title("standard deviations of histograms, df = %g" % df)
    plt.show()
    stop = time.time()
    print(stop-start)
    

sigma_plotter_mat()