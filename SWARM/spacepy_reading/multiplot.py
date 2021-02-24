
import SWARMprocess
import numpy as np
import matplotlib.pyplot as plt
import time
pro = SWARMprocess.SWARMprocess()

def sigma_plotter():
    f0 = 0.1
    f1 = 0.9
    n = 400
    df = 0.05
    start = time.time()
    object = multi_SWARMreader.MultiSWARM(2013, 12, 9, 2013, 12, 31)
    freq0s, sigmas, means = object.freq_sig(df = df, n = n, f0 = f0, f1 = f1, bins_ = 50)

    plt.plot(freq0s, sigmas[0])
    plt.plot(freq0s, sigmas[1])
    plt.plot(freq0s, sigmas[2])
    plt.xlabel("Lower integral limit")
    plt.ylabel("STD")
    plt.legend(["B-A", "B-C", "A-C"])
    plt.title("standard deviations of histograms, df = 0.01")
    plt.show()
    stop = time.time()
    print(stop-start)


def sigma_plotter_lat():
    df = 0.05
    f0 = 0.1
    f1 = 0.9
    pole = "north"
    start = time.time()
    object = multi_SWARMreader.MultiSWARM(2013, 12, 9, 2013, 12, 25)
    freq0s, sigmas, means = object.freq_sig_lat(df = df, n = 400, bins_ = 50, lat1 = 90, lat0 = 75,
                                                f0 = f0, f1 = f1, pole = pole)

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

def spec_lat_hist_plotter():
    start = time.time()
    object = multi_SWARMreader.MultiSWARM(2013, 12, 14, 2013, 12, 14)
    minfreq = 0.1
    maxfreq = 0.3
    hists, bins = object.multi_histmake_lat(bins_ = 50, minfreq = minfreq,\
                                            maxfreq = maxfreq, lat0 = 65, lat1 = 77,\
                                            norm = False, n = 300, mlat = True)



    means = np.zeros(len(bins))
    stds = np.zeros_like(means)
    for i in range(len(bins)):
        dbins = bins[i][1]-bins[i][0]
        hists[i] = hists[i]/np.sum(hists[i]*dbins)
        curr_std, curr_mean = object.pro.std_mean(hists[i], bins[i])
        means[i] = curr_mean
        stds[i] = curr_std


    xs = np.linspace(np.min(bins[0]), np.max(bins[0]), 1000)

    fig, axs = plt.subplots(2, 3)

    for j in range(3):
        for i in [0, 1]:
            k = 2*j+i
            gauss = object.pro.gauss_curve(x = xs, mean = means[k], std = stds[k])
            axs[i, j].plot(xs, gauss, "r")
            width = bins[k][1] - bins[k][0]
            axs[i, j].bar(bins[k], hists[k], width = 0.9*width)
            # axs[i, j].xlabel("Relative difference")
            # axs[i, j].ylabel("Normalized occurence")
            # axs[i, j].legend(["Gauss approximation", "Histogram"])
            stringies = ["B-A high", "B-A low", "B-C high", "B-C low", "A-C high", "A-C low"]
            axs[i, j].set_title("std = %g" % stds[k] +", " + stringies[k])
    for ax in axs.flat:
        ax.set(xlabel="relative diff", ylabel="normalized occurence")

    # for i in range(len(hists)):
    #     gauss = object.pro.gauss_curve(x = xs, mean = means[i], std = stds[i])
    #     plt.figure(i)
    #     width = bins[i][1] - bins[i][0]
    #     plt.plot(xs, gauss, "r")
    #     plt.bar(bins[i], hists[i], width = 0.9*width)
    #     plt.xlabel("Relative difference")
    #     plt.ylabel("Normalized occurence")
    #     plt.legend(["Gauss approximation", "Histogram"])
    #     stringies = ["B-A high", "B-A low", "B-C high", "B-C low", "A-C high", "A-C low"]
    #     plt.title(stringies[i])

    plt.show()

    # BA_shift = object.BA_shift
    # BC_shift = object.BC_shift
    # AC_shift = BC_shift - BA_shift
    # BA_shift = BA_shift/2
    # BC_shift = BC_shift/2
    # AC_shift = AC_shift/2
    # plt.plot([BA_shift, AC_shift, BC_shift], [stds[0], stds[4], stds[2]], "o-")
    # plt.xlabel("Time difference between satellites")
    # plt.ylabel("standard deviation in high latitude histograms")
    # plt.show()


    stop = time.time()
    print(stop-start)

def std_timeshift():
    start = time.time()
    minfreq = 0.1
    maxfreq = 0.3
    day0 = 9
    day1 = 31
    lat1 = 77
    lat0 = 65
    mlat = True
    shift_list = []
    std_list = []
    k_high = []
    for day in range(day0, day1):
        object = multi_SWARMreader.MultiSWARM(2013, 12, day, 2013, 12, day)
        hists, bins = object.multi_histmake_lat(bins_ = 50, minfreq = minfreq,\
                                                maxfreq = maxfreq, lat0 = lat0, lat1 = lat1,\
                                                norm = False, n = 300, pole = "north", mlat = mlat)
        if object.samelength != True:
            continue



        means = np.zeros(len(bins))
        stds = np.zeros_like(means)
        for i in range(len(bins)):
            dbins = bins[i][1]-bins[i][0]
            hists[i] = hists[i]/np.sum(hists[i]*dbins)
            curr_std, curr_mean = object.pro.std_mean(hists[i], bins[i])
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
        std_list.append(stds[2])
        shift_list.append(AC_shift)
        std_list.append(stds[4])
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

def std_timeshift_n():
    start = time.time()
    minfreq = 0.5
    maxfreq = 0.7
    day0 = 9
    day1 = 31
    lat1 = 90
    lat0 = 75
    mlat = True


    a_list_list = []
    b_list_list = []

    minfreqs = [0.1, 0.3, 0.5, 0.7]
    maxfreqs = [0.3, 0.5, 0.7, 0.9]


    nlist = np.linspace(10, 600, 10)
    for f in range(len(minfreqs)):
        a_list = []
        b_list = []
        for n in nlist:
            print("work is %g%% done" % ((n/nlist[-1]/len(minfreqs) + f/(len(minfreqs)))*100))
            shift_list = []
            std_list = []
            k_high = []
            for day in range(day0, day1):
                object = multi_SWARMreader.MultiSWARM(2013, 12, day, 2013, 12, day)
                hists, bins = object.multi_histmake_lat(bins_ = 50, minfreq = minfreqs[f],\
                                                        maxfreq = maxfreqs[f], lat0 = lat0, lat1 = lat1,\
                                                        norm = True, n = n, pole = "north", mlat = mlat)
                if object.samelength != True:
                    continue



                means = np.zeros(len(bins))
                stds = np.zeros_like(means)
                for i in range(len(bins)):
                    dbins = bins[i][1]-bins[i][0]
                    hists[i] = hists[i]/np.sum(hists[i]*dbins)
                    curr_std, curr_mean = object.pro.std_mean(hists[i], bins[i])
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
                std_list.append(stds[2])
                shift_list.append(AC_shift)
                std_list.append(stds[4])
                if day == 14 or day == 25:
                    for j in range(3):
                        k_high.append(1)
                else:
                    for j in range(3):
                        k_high.append(0)

            p = np.polyfit(shift_list, std_list, deg = 1)
            a = p[0]; b = p[1]
            a_list.append(a)
            b_list.append(b)
        a_list = np.array(a_list)
        b_list = np.array(b_list)
        a_list_list.append(a_list)
        b_list_list.append(b_list)


    legend = []
    plt.figure(0)
    for i in range(len(a_list_list)):
        plt.plot(nlist/2, a_list_list[i])
        legend.append("f: %g - %g" % (minfreqs[i], maxfreqs[i]))

    plt.xlabel("Window size [s]")
    plt.ylabel("linear regression coefficient a")
    plt.title("temporal regression component as function of window size")
    plt.legend(legend)

    plt.figure(1)
    for i in range(len(b_list_list)):
        plt.plot(nlist/2, b_list_list[i])

    plt.xlabel("Window size [s]")
    plt.ylabel("linear regression constant b")
    plt.title("spatial regression component as function of window size")
    plt.legend(legend)
    plt.show()




    stop = time.time()
    print(stop-start)


def multi_std_timeshift():
    start = time.time()
    minfreq = 0.1
    maxfreq = 0.3
    day0 = 9
    day1 = 31
    lat1 = 90
    lat0 = 77
    mlat = True
    ns = [80, 300]
    colors = ["k", "r"]
    for k in range(len(ns)):
        shift_list = []
        std_list = []
        k_high = []
        for day in range(day0, day1):
            object = multi_SWARMreader.MultiSWARM(2013, 12, day, 2013, 12, day)
            hists, bins = object.multi_histmake_lat(bins_ = 50, minfreq = minfreq,\
                                                    maxfreq = maxfreq, lat0 = lat0, lat1 = lat1,\
                                                    norm = True, n = ns[k], pole = "north", mlat = mlat)
            if object.samelength != True:
                continue



            means = np.zeros(len(bins))
            stds = np.zeros_like(means)
            for i in range(len(bins)):
                dbins = bins[i][1]-bins[i][0]
                hists[i] = hists[i]/np.sum(hists[i]*dbins)
                curr_std, curr_mean = object.pro.std_mean(hists[i], bins[i])
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
            std_list.append(stds[2])
            shift_list.append(AC_shift)
            std_list.append(stds[4])
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
        xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
        plt.figure(1)
        plt.plot(xs, xs*a + b, colors[k])
        plt.plot(shift_list, std_list, colors[k] + "o")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("Standard deviation of histograms")
    plt.title("Stds per day, f: %g - %g" % (minfreq, maxfreq))
    plt.legend(["n = 80", "n = 80", "n = 300", "n = 300"])
    plt.show()

    stop = time.time()
    print(stop-start)
