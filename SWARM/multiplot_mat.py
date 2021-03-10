import multi_matreader
import general_matreader
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
    lat0 = 7
    n = 100
    shift_list = []
    std_list = []
    k_high = []

    for day in range(day0, day1+1):
        object = multi_matreader.MultiMat(day, day)
        hists, bins = object.multi_histmake(n = n, minfreq = minfreq,\
        maxfreq = maxfreq, bins_ = 200, lat1 = lat1, lat0 = lat0,\
        norm = True, pole = "both")
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
    df = 0.05
    jump = 0.01
    lat1 = 90
    lat0 = 77
    start = time.time()
    object = multi_matreader.MultiMat(9, 31)
    freq0s, sigmas, means = object.freq_sig(df = df, jump = jump, n = n,\
                                            f0 = f0, f1 = f1, bins_ = 50,\
                                            abs = False, norm = False,\
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


def plothing():
    """
    random plot
    """
    multi = multi_matreader.MultiMat(21, 21)
    object = general_matreader.MatReader(multi.filepath(multi.init_loop_index))
    start = 1325*2 #index of 16 minutes
    stop = 1375*2 #index of 26 minutes
    seconds = object.secondsA[start:stop]
    NeA = object.NeA[start:stop]
    plt.figure(0)
    plt.plot(seconds, NeA)
    plt.xlabel("Time [s]")
    plt.ylabel("NeA")
    plt.title("Raw data from sat A")
    # plt.show()

    plt.figure(1)
    NeA = pro.meanie(NeA, 5)
    plt.plot(seconds, NeA)
    plt.xlabel("Time [s]")
    plt.ylabel("Ne")
    plt.title("Data, rolling mean")
    # plt.show()

    fs = 2
    N = len(NeA)
    freqs = np.linspace(-fs/2, fs/2, N)
    fft = np.fft.fft(NeA)/N
    fft = np.roll(fft, int(N/2))
    plt.figure(2)
    plt.plot(freqs, np.log10(np.abs(fft)))
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Log10(Fourier Coefficients)")
    plt.title("Power Spectral Density")
    plt.show()

def histplot():
    """
    plots histograms
    """
    day0 = 9
    day1 = 31
    object = multi_matreader.MultiMat(day0, day1)
    n = 100
    minfreq = 0.1
    maxfreq = 0.3
    bins_ = 100
    lat1 = 90
    lat0 = 77
    hists, bins = object.multi_histmake(n, minfreq, maxfreq, bins_, lat1, lat0,\
                       abs = False, norm = True, pole = "north")

    stds = np.zeros(len(hists))
    means = np.zeros_like(stds)
    width = bins[1] - bins[0]
    for j in range(len(hists)):
        hists[j] = hists[j]/np.sum(hists[j]*width)
        curr_std, curr_mean = object.pro.std_mean(hists[j], bins)
        stds[j] = curr_std
        means[j] = curr_mean

    xs = np.linspace(-1, 1, 1000)
    gauss1 = object.pro.gauss_curve(xs, means[0], stds[0])
    gauss2 = object.pro.gauss_curve(xs, means[1], stds[1])
    gauss3 = object.pro.gauss_curve(xs, means[2], stds[2])

    fig, axs = plt.subplots(1, 3)
    axs[0].set_title("B - A, $\sigma$ = %g, $\mu$ = %g" % (stds[0], means[0]))
    axs[0].plot(xs, gauss1, "r")
    axs[0].bar(bins, hists[0], width = width)

    axs[1].set_title("B - C, $\sigma$ = %g, $\mu$ = %g" % (stds[1], means[1]))
    axs[1].plot(xs, gauss2, "r")
    axs[1].bar(bins, hists[1], width = width)

    axs[2].set_title("A - C, $\sigma$ = %g, $\mu$ = %g" % (stds[2], means[2]))
    axs[2].plot(xs, gauss3, "r")
    axs[2].bar(bins, hists[2], width = width)


    xlabels = ["relative difference","relative difference","relative difference"]
    ylabels = ["Norm occ","Norm occ","Norm occ"]
    for i in range(len(axs.flat)):
        axs.flat[i].set(xlabel=xlabels[i], ylabel=ylabels[i])
        #axs.flat[i].set_aspect("equal", "box")
        #axs.flat[i].set_xlim(-1, 1)
        #axs.flat[i].set_ylim(0, 5)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()

    plt.legend(["normal distribution", "data"])
    plt.show()

std_timeshift_mat()