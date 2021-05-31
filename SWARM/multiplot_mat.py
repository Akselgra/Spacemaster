import multi_matreader
import general_matreader
import SWARMprocess
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.stats import skew
from scipy.stats import kurtosis
pro = SWARMprocess.SWARMprocess()

fig_width_pt = 418.0  # Get this from LaTeX using \showthe\columnwidth
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


def std_timeshift_mat():
    start = time.time()
    minfreq = 0.1
    maxfreq = 1
    day0 = 9
    day1 = 31
    lat1 = -90
    lat0 = -77
    n = 120
    shift_list = []
    std_list = []
    k_high = []

    for day in range(day0, day1+1):
        object = multi_matreader.MultiMat(day, day)
        hists, bins = object.multi_histmake(n = n, minfreq = minfreq,\
        maxfreq = maxfreq, bins_ = 200, lat1 = lat1, lat0 = lat0,\
        norm = True, pole = "south")
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



    # p = np.polyfit(shift_list, std_list, deg = 1)
    # a = p[0]; b = p[1]
    a, b, adiv, bdiv = pro.linear_regression(shift_list, std_list)
    print("Slope of regression is %g pm %g" % (a, adiv))
    print("Constant of linear regression is %g pm %g" % (b, bdiv))
    print("a is %g standard deviations away from zero" % (a/adiv))

    print(a/b*100)
    xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
    plt.figure(1)
    plt.plot(xs, xs*a + b)
    plt.plot(shift_list, std_list, "ko")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("Standard deviation of histograms")
    plt.title("Stds per day, f: %g - %g, a = %g, b = %g" % (minfreq, maxfreq, a, b))
    plt.legend(["Linear regression", "data points"])
    plt.savefig("Figures/matfigs/std_timeshift/std_timeshift_polarcap_south.pdf")
    
    plt.show()
    
    
    # plt.figure(2)
    # for i in range(len(shift_list)):
    #     if k_high[i] == 0:
    #         plt.plot(shift_list[i], std_list[i], "ko")
    #     else:
    #         plt.plot(shift_list[i], std_list[i], "ro")


    # plt.xlabel("Time between satellites [s]")
    # plt.ylabel("Standard deviation of histograms")
    # plt.title("Stds per day, integral limits: %g to %g" % (minfreq, maxfreq))
    # plt.show()



    stop = time.time()
    print(stop-start)

def skewness_timeshift_mat():
    start = time.time()
    minfreq = 0.1
    maxfreq = 1
    day0 = 9
    day1 = 31
    lat1 = 90
    lat0 = 77
    n = 120
    shift_list = []
    skew_list = []
    kurtosis_list = []
    mean_list = []
    std_list = []

    region = "Northpole"
    for day in range(day0, day1+1):
        object = multi_matreader.MultiMat(day, day)
        
        I_BA, I_BC, I_AC = object.multi_compind(n = n, minfreq = minfreq,\
        maxfreq = maxfreq, lat1 = lat1, lat0 = lat0,\
        norm = True, pole = "north")
        
        # AC_std = np.std(I_AC)
        # AC_mean = np.mean(I_AC)
        # xs = np.linspace(-1, 1, 1000)
        # gauss = object.pro.gauss_curve(xs, AC_mean, AC_std)
        # AChist, bins = np.histogram(I_AC)
        # bins = (bins[1:] + bins[:-1])/2
        # width = bins[1] - bins[0]
        # AChist = AChist/np.sum(AChist*width)
        # plt.bar(bins, AChist, width = width)
        # plt.plot(xs, gauss, "r")
        # plt.show()
        
        curr_skew_BA = skew(I_BA)
        curr_skew_BC = skew(I_BC)
        curr_skew_AC = skew(I_AC)
        
        curr_kurtosis_BA = kurtosis(I_BA)
        curr_kurtosis_BC = kurtosis(I_BC)
        curr_kurtosis_AC = kurtosis(I_AC)
        
        curr_mean_BA = np.mean(I_BA)
        curr_mean_BC = np.mean(I_BC)
        curr_mean_AC = np.mean(I_AC)
        
        curr_std_BA = np.std(I_BA)
        curr_std_BC = np.std(I_BC)
        curr_std_AC = np.std(I_AC)
        
        
        
        

        BA_shift = object.BA_shift
        BC_shift = object.BC_shift
        AC_shift = BC_shift - BA_shift
        BA_shift = BA_shift/2
        BC_shift = BC_shift/2
        AC_shift = AC_shift/2
        
        shift_list.append(BA_shift)
        skew_list.append(curr_skew_BA)
        kurtosis_list.append(curr_kurtosis_BA)
        mean_list.append(curr_mean_BA)
        std_list.append(curr_std_BA)
        
        shift_list.append(BC_shift)
        skew_list.append(curr_skew_BC)
        kurtosis_list.append(curr_kurtosis_BC)
        mean_list.append(curr_mean_BC)
        std_list.append(curr_std_BC)
        
        shift_list.append(AC_shift)
        skew_list.append(curr_skew_AC)
        kurtosis_list.append(curr_kurtosis_AC)
        mean_list.append(curr_mean_AC)
        std_list.append(curr_std_AC)



    # p = np.polyfit(shift_list, std_list, deg = 1)
    # a = p[0]; b = p[1]
    a, b, adiv, bdiv = pro.linear_regression(shift_list, skew_list)
    print("Slope of regression is %g pm %g" % (a, adiv))
    print("Constant of linear regression is %g pm %g" % (b, bdiv))
    print("a is %g standard deviations away from zero" % (a/adiv))

    print(a/b*100)
    xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
    plt.figure(1)
    plt.ylim([-3.2, 3.2])
    plt.plot(xs, xs*a + b)
    plt.plot(shift_list, skew_list, "ko")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("Skewness of histograms")
    plt.title("Skewness per day, f: %g - %g, a = %g, b = %g" % (minfreq, maxfreq, a, b))
    plt.legend(["Linear regression", "data points"])
    plt.savefig("Figures/matfigs/skew/skew_" + region + ".pdf")
    plt.show()
    
    
    # p = np.polyfit(shift_list, std_list, deg = 1)
    # a = p[0]; b = p[1]
    a, b, adiv, bdiv = pro.linear_regression(shift_list, kurtosis_list)
    print("Slope of regression is %g pm %g" % (a, adiv))
    print("Constant of linear regression is %g pm %g" % (b, bdiv))
    print("a is %g standard deviations away from zero" % (a/adiv))

    print(a/b*100)
    xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
    plt.figure(1)
    plt.ylim([-5, 22])
    plt.plot(xs, xs*a + b)
    plt.plot(shift_list, kurtosis_list, "ko")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("kurtosis of histograms")
    plt.title("Kurtosis per day, f: %g - %g, a = %g, b = %g" % (minfreq, maxfreq, a, b))
    plt.legend(["Linear regression", "data points"])
    plt.savefig("Figures/matfigs/kurtosis/kurtosis_" + region + ".pdf")
    plt.show()
    
    a, b, adiv, bdiv = pro.linear_regression(shift_list, mean_list)
    print("Slope of regression is %g pm %g" % (a, adiv))
    print("Constant of linear regression is %g pm %g" % (b, bdiv))
    print("a is %g standard deviations away from zero" % (a/adiv))

    print(a/b*100)
    xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
    plt.figure(1)
    plt.ylim([-0.15, 0.15])
    plt.plot(xs, xs*a + b)
    plt.plot(shift_list,  mean_list, "ko")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("Mean of histograms")
    plt.title("Mean per day, f: %g - %g, a = %g, b = %g" % (minfreq, maxfreq, a, b))
    plt.legend(["Linear regression", "data points"])
    plt.savefig("Figures/matfigs/mean/mean_" + region + ".pdf")
    plt.show()
    
    
    a, b, adiv, bdiv = pro.linear_regression(shift_list, std_list)
    print("Slope of regression is %g pm %g" % (a, adiv))
    print("Constant of linear regression is %g pm %g" % (b, bdiv))
    print("a is %g standard deviations away from zero" % (a/adiv))

    print(a/b*100)
    xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
    plt.figure(1)
    plt.ylim([0, 0.35])
    plt.plot(xs, xs*a + b)
    plt.plot(shift_list,  std_list, "ko")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("Standard deviation of histograms")
    plt.title("std per day, f: %g - %g, a = %g, b = %g" % (minfreq, maxfreq, a, b))
    plt.legend(["Linear regression", "data points"])
    plt.savefig("Figures/matfigs/std/std_" + region + ".pdf")
    plt.show()



    stop = time.time()
    print(stop-start)
def sigma_plotter_mat():
    """
    plots std as a function of frequency
    """
    f0 = 0
    f1 = 1
    n = 120
    df = 0.5
    jump = 0.01
    lat1 = 90
    lat0 = 77
    start = time.time()
    object = multi_matreader.MultiMat(15, 15)
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
    plt.grid("on")
    plt.savefig("Figures/matfigs/sigma_f/mat_sigma_f.pdf")
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

def multiplot_moments():
    start = time.time()
    minfreq = 0.1
    maxfreq = 1
    day0 = 9
    day1 = 31
    lat1 = 30
    lat0 = 0
    n = 120
    

    region = "equatorial"
    lat0s = [0, 30, 70, 77]
    lat1s = [30, 70, 77, 90]
    regions = ["Equatorial (A)", "Midlatitude (B)", "Auroral Oval (C)", "Northern polar cap (D)"]
    poles = ["both", "both", "north", "north"]
    shift_list_list = []
    kurtosis_list_list = []
    skew_list_list = []
    std_list_list = []
    mean_list_list = []
    
    for l in range(len(lat0s)):
        shift_list = []
        skew_list = []
        kurtosis_list = []
        mean_list = []
        std_list = []
        for day in range(day0, day1+1):
            object = multi_matreader.MultiMat(day, day)
            
            I_BA, I_BC, I_AC = object.multi_compind(n = n, minfreq = minfreq,\
            maxfreq = maxfreq, lat1 = lat1s[l], lat0 = lat0s[l],\
            norm = True, pole = poles[l])
            

            
            curr_skew_BA = skew(I_BA)
            curr_skew_BC = skew(I_BC)
            curr_skew_AC = skew(I_AC)
            
            curr_kurtosis_BA = kurtosis(I_BA)
            curr_kurtosis_BC = kurtosis(I_BC)
            curr_kurtosis_AC = kurtosis(I_AC)
            
            curr_mean_BA = np.mean(I_BA)
            curr_mean_BC = np.mean(I_BC)
            curr_mean_AC = np.mean(I_AC)
            
            curr_std_BA = np.std(I_BA)
            curr_std_BC = np.std(I_BC)
            curr_std_AC = np.std(I_AC)
            
            
            
            
    
            BA_shift = object.BA_shift
            BC_shift = object.BC_shift
            AC_shift = BC_shift - BA_shift
            BA_shift = BA_shift/2
            BC_shift = BC_shift/2
            AC_shift = AC_shift/2
            
            shift_list.append(BA_shift)
            skew_list.append(curr_skew_BA)
            kurtosis_list.append(curr_kurtosis_BA)
            mean_list.append(curr_mean_BA)
            std_list.append(curr_std_BA)
            
            shift_list.append(BC_shift)
            skew_list.append(curr_skew_BC)
            kurtosis_list.append(curr_kurtosis_BC)
            mean_list.append(curr_mean_BC)
            std_list.append(curr_std_BC)
            
            shift_list.append(AC_shift)
            skew_list.append(curr_skew_AC)
            kurtosis_list.append(curr_kurtosis_AC)
            mean_list.append(curr_mean_AC)
            std_list.append(curr_std_AC)
    
    
    
        shift_list_list.append(shift_list)
        kurtosis_list_list.append(kurtosis_list)
        skew_list_list.append(skew_list)
        std_list_list.append(std_list)
        mean_list_list.append(mean_list)
        
    
    fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
    
    print("\\begin{table}[htbp]")
    print("\\centering")
    print("\\caption{Table of regression coefficients for the linear regressions found in Figure \\ref{fig:multi_kurtosis}. The first column has the latitudinal regions, the second column has the regression slopes and the third column has the regression constants.}")
    print("\\begin{tabular}{|c|c|c|}")
    print("\\hline")
    print("Region & Kurtosis regression slope [s$^{-1}$] & Kurtosis regression constant [s$^{-1}$]\\" + "\\")
    print("\\hline")
    for i in range(2):
        for j in range(2):
            l = (i*2 + j)
            a, b, adiv, bdiv = pro.linear_regression(shift_list_list[l], kurtosis_list_list[l])
            a = float(a); b = float(b); adiv = float(adiv); bdiv = float(bdiv)
            
            norm_a = pro.standard_form(a, dec = 1)
            norm_b = pro.standard_form(b, dec = 1)
            norm_adiv = pro.standard_form(adiv, dec = 1)
            norm_bdiv = pro.standard_form(bdiv, dec = 1)
            print(regions[l] + " & $" + norm_a + " \pm " + norm_adiv + "$ & $" + norm_b + " \pm " + norm_bdiv + "$\\" + "\\")
            print("\\hline")
            
            # print("region = " + regions[l])
            # print("Slope of regression is %g pm %g" % (a, adiv))
            # print("Constant of linear regression is %g pm %g" % (b, bdiv))
            xs = np.linspace(np.min(shift_list_list[l]), np.max(shift_list_list[l]), 1000)
            axs[i, j].plot(shift_list_list[l], kurtosis_list_list[l], "k.")
            axs[i, j].plot(xs, xs*a + b)
            axs[i, j].set_title(regions[l])
            axs[i, j].grid("on")
    
    for ax in axs.flat:
        ax.set(xlabel = "Time difference $\Delta t$ [s]", ylabel = "Kurtosis")
    
    for ax in axs.flat:
        ax.label_outer()
    
    print("\\end{tabular}")
    print("\\label{tab:multi_kurtosis}")
    print("\\end{table}")
    plt.savefig("Figures/matfigs/kurtosis/multi_kurtosis.pdf")
    plt.show()
    
    #-----------------------------------------------------------
    
    print("\\begin{table}[htbp]")
    print("\\centering")
    print("\\caption{Table of regression coefficients for the linear regressions found in Figure \\ref{fig:multi_skewness}. The first column has the latitudinal regions, the second column has the regression slopes and the third column has the regression constants.}")
    print("\\begin{tabular}{|c|c|c|}")
    print("\\hline")
    print("Region & Skewness regression slope [s$^{-1}$] & Skewness regression constant [s$^{-1}$]\\" + "\\")
    print("\\hline")
    
    fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
    for i in range(2):
        for j in range(2):
            l = (i*2 + j)
            a, b, adiv, bdiv = pro.linear_regression(shift_list_list[l], skew_list_list[l])
            a = float(a); b = float(b); adiv = float(adiv); bdiv = float(bdiv)
            
            norm_a = pro.standard_form(a, dec = 1)
            norm_b = pro.standard_form(b, dec = 1)
            norm_adiv = pro.standard_form(adiv, dec = 1)
            norm_bdiv = pro.standard_form(bdiv, dec = 1)
            print(regions[l] + " & $" + norm_a + " \pm " + norm_adiv + "$ & $" + norm_b + " \pm " + norm_bdiv + "$\\" + "\\")
            print("\\hline")
            
            # print("region = " + regions[l])
            # print("Slope of regression is %g pm %g" % (a, adiv))
            # print("Constant of linear regression is %g pm %g" % (b, bdiv))
            xs = np.linspace(np.min(shift_list_list[l]), np.max(shift_list_list[l]), 1000)
            axs[i, j].plot(shift_list_list[l], skew_list_list[l], "k.")
            axs[i, j].plot(xs, xs*a + b)
            axs[i, j].set_title(regions[l])
            axs[i, j].grid("on")
    
    for ax in axs.flat:
        ax.set(xlabel = "Time difference $\Delta t$ [s]", ylabel = "Skewness")
    
    for ax in axs.flat:
        ax.label_outer()
    
    print("\\end{tabular}")
    print("\\label{tab:multi_skewness}")
    print("\\end{table}")
    plt.savefig("Figures/matfigs/skew/multi_skewness.pdf")
    plt.show()
    
    #---------------------------------------------------------------------
    
    print("\\begin{table}[htbp]")
    print("\\centering")
    print("\\caption{Table of regression coefficients for the linear regressions found in Figure \\ref{fig:multi_std}. The first column has the latitudinal regions, the second column has the regression slopes and the third column has the regression constants.}")
    print("\\begin{tabular}{|c|c|c|}")
    print("\\hline")
    print("Region & Std regression slope [s$^{-1}$] & Std regression constant [s$^{-1}$]\\" + "\\")
    print("\\hline")
    
    fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
    for i in range(2):
        for j in range(2):
            l = (i*2 + j)
            a, b, adiv, bdiv = pro.linear_regression(shift_list_list[l], std_list_list[l])
            a = float(a); b = float(b); adiv = float(adiv); bdiv = float(bdiv)
            
            norm_a = pro.standard_form(a, dec = 1)
            norm_b = pro.standard_form(b, dec = 1)
            norm_adiv = pro.standard_form(adiv, dec = 1)
            norm_bdiv = pro.standard_form(bdiv, dec = 1)
            print(regions[l] + " & $" + norm_a + " \pm " + norm_adiv + "$ & $" + norm_b + " \pm " + norm_bdiv + "$\\" + "\\")
            print("\\hline")
            # print("region = " + regions[l])
            # print("Slope of regression is %g pm %g" % (a, adiv))
            # print("Constant of linear regression is %g pm %g" % (b, bdiv))
            xs = np.linspace(np.min(shift_list_list[l]), np.max(shift_list_list[l]), 1000)
            axs[i, j].plot(shift_list_list[l], std_list_list[l], "k.")
            axs[i, j].plot(xs, xs*a + b)
            axs[i, j].set_title(regions[l])
            axs[i, j].grid("on")
    
    for ax in axs.flat:
        ax.set(xlabel = "Time difference $\Delta t$ [s]", ylabel = "Standard deviation")
    
    for ax in axs.flat:
        ax.label_outer()
        
    print("\\end{tabular}")
    print("\\label{tab:multi_std}")
    print("\\end{table}")
    
    plt.savefig("Figures/matfigs/std/multi_std.pdf")
    plt.show()
    
    
    #-----------------------------------------------------------------------
    
    
    print("\\begin{table}[htbp]")
    print("\\centering")
    print("\\caption{Table of regression coefficients for the linear regressions found in Figure \\ref{fig:multi_mean}. The first column has the latitudinal regions, the second column has the regression slopes and the third column has the regression constants.}")
    print("\\begin{tabular}{|c|c|c|}")
    print("\\hline")
    print("Region & Mean regression slope [s$^{-1}$] & Mean regression constant [s$^{-1}$]\\" + "\\")
    print("\\hline")
    
    fig, axs = plt.subplots(2, 2, sharex = True, sharey = True)
    for i in range(2):
        for j in range(2):
            l = (i*2 + j)
            a, b, adiv, bdiv = pro.linear_regression(shift_list_list[l], mean_list_list[l])
            a = float(a); b = float(b); adiv = float(adiv); bdiv = float(bdiv)
            
            norm_a = pro.standard_form(a, dec = 1)
            norm_b = pro.standard_form(b, dec = 1)
            norm_adiv = pro.standard_form(adiv, dec = 1)
            norm_bdiv = pro.standard_form(bdiv, dec = 1)
            print(regions[l] + " & $" + norm_a + " \pm " + norm_adiv + "$ & $" + norm_b + " \pm " + norm_bdiv + "$\\" + "\\")
            print("\\hline")
            # print("region = " + regions[l])
            # print("Slope of regression is %g pm %g" % (a, adiv))
            # print("Constant of linear regression is %g pm %g" % (b, bdiv))
            xs = np.linspace(np.min(shift_list_list[l]), np.max(shift_list_list[l]), 1000)
            axs[i, j].plot(shift_list_list[l], mean_list_list[l], "k.")
            axs[i, j].plot(xs, xs*a + b)
            axs[i, j].set_title(regions[l])
            axs[i, j].grid("on")
    
    for ax in axs.flat:
        ax.set(xlabel = "Time difference $\Delta t$ [s]", ylabel = "Mean")
    
    for ax in axs.flat:
        ax.label_outer()
    
    print("\\end{tabular}")
    print("\\label{tab:multi_mean}")
    print("\\end{table}")
    plt.savefig("Figures/matfigs/mean/multi_mean.pdf")
    plt.show()
    
    
        
    
    stop = time.time()
    print(stop-start)
def histplot():
    """
    plots histograms
    """
    day0 = 15
    day1 = 15
    object = multi_matreader.MultiMat(day0, day1)
    n = 120
    minfreq = 0.1
    maxfreq = 1
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
    print(means[0])
    print(stds[0])
    print(object.BA_shift)
    gauss1 = object.pro.gauss_curve(xs, means[0], stds[0])
    gauss2 = object.pro.gauss_curve(xs, means[1], stds[1])
    gauss3 = object.pro.gauss_curve(xs, means[2], stds[2])

    # fig, axs = plt.subplots(1, 3)
    # axs[0].set_title("B - A, $\sigma$ = %g, $\mu$ = %g" % (stds[0], means[0]))
    # axs[0].plot(xs, gauss1, "r")
    # axs[0].bar(bins, hists[0], width = width)

    # axs[1].set_title("B - C, $\sigma$ = %g, $\mu$ = %g" % (stds[1], means[1]))
    # axs[1].plot(xs, gauss2, "r")
    # axs[1].bar(bins, hists[1], width = width)

    # axs[2].set_title("A - C, $\sigma$ = %g, $\mu$ = %g" % (stds[2], means[2]))
    # axs[2].plot(xs, gauss3, "r")
    # axs[2].bar(bins, hists[2], width = width)


    # xlabels = ["relative difference","relative difference","relative difference"]
    # ylabels = ["Norm occ","Norm occ","Norm occ"]
    # for i in range(len(axs.flat)):
    #     axs.flat[i].set(xlabel=xlabels[i], ylabel=ylabels[i])
    #     #axs.flat[i].set_aspect("equal", "box")
    #     #axs.flat[i].set_xlim(-1, 1)
    #     #axs.flat[i].set_ylim(0, 5)

    # # Hide x labels and tick labels for top plots and y ticks for right plots.
    # for ax in axs.flat:
    #     ax.label_outer()

    # plt.legend(["normal distribution", "data"])
    # plt.show()
    
    plt.bar(bins, hists[0], width = width)
    plt.plot(xs, gauss1, "r")
    plt.xlabel("Comparison index $I_{B-A}$")
    plt.ylabel("Normalized occurence")
    plt.legend(["Gaussian distribution", "Data"])
    plt.title("Histogram of comparison indices for polar cap region")
    plt.savefig("Figures/matfigs/histograms/I_BA_pole_hist.pdf")
    plt.show()

def std_distance():
    start = time.time()
    minfreq = 0.1
    maxfreq = 0.9
    day0 = 9
    day1 = 31
    lat1 = 90
    lat0 = 1
    n = 150
    lat_diff_list = []
    long_diff_list = []
    dist_list_lat = []
    dist_list_long = []
    dist_list = []
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

        latA = object.latA
        latB = object.latB
        latC = object.latC

        longA = object.longA
        longB = object.longB
        longC = object.longC

        radA = object.radA
        radB = object.radB
        radC = object.radC

        lat_dist_BA = np.deg2rad(np.abs(np.mean(latB - latA)))*np.mean(object.radB)
        lat_dist_BC = np.deg2rad(np.abs(np.mean(latB - latC)))*np.mean(object.radB)
        lat_dist_AC = np.deg2rad(np.abs(np.mean(latA - latC)))*np.mean(object.radA)

        long_dist_BA = np.deg2rad(np.abs(np.mean(longB - longA)))*np.mean(object.radB)
        long_dist_BC = np.deg2rad(np.abs(np.mean(longB - longC)))*np.mean(object.radB)
        long_dist_AC = np.deg2rad(np.abs(np.mean(longA - longC)))*np.mean(object.radA)

        BA_dist = np.mean(pro.great_circle_distance(latB, longB, radB, latA, longA, radA))
        BC_dist = np.mean(pro.great_circle_distance(latB, longB, radB, latC, longC, radC))
        AC_dist = np.mean(pro.great_circle_distance(latA, longA, radA, latC, longC, radC))


        BA_shift = BA_shift/2
        BC_shift = BC_shift/2
        AC_shift = AC_shift/2
        shift_list.append(BA_shift)
        std_list.append(stds[0])
        shift_list.append(BC_shift)
        std_list.append(stds[1])
        shift_list.append(AC_shift)
        std_list.append(stds[2])

        lat_diff_list.append(np.mean(np.abs(latB - latA)))
        lat_diff_list.append(np.mean(np.abs(latB - latC)))
        lat_diff_list.append(np.mean(np.abs(latA - latC)))
        
        

        long_diff_list.append(np.mean(np.abs(longB - longA)))
        long_diff_list.append(np.mean(np.abs(longB - longC)))
        long_diff_list.append(np.mean(np.abs(longA - longC)))
        
        

        dist_list_lat.append(lat_dist_BA)
        dist_list_lat.append(lat_dist_BC)
        dist_list_lat.append(lat_dist_AC)

        dist_list_long.append(long_dist_BA)
        dist_list_long.append(long_dist_BC)
        dist_list_long.append(long_dist_AC)

        dist_list.append(BA_dist)
        dist_list.append(BC_dist)
        dist_list.append(AC_dist)


        if day == 14 or day == 25:
            for j in range(3):
                k_high.append(1)
        else:
            for j in range(3):
                k_high.append(0)


    # p = np.polyfit(shift_list, std_list, deg = 1)
    # a = p[0]; b = p[1]
    # print("Slope of regression is %g" % a)
    # print("Constant of linear regression is %g" % b)

    # print(a/b*100)
    # xs = np.linspace(np.min(shift_list), np.max(shift_list), 1000)
    plt.figure(1)
    # plt.plot(xs, xs*a + b)
    plt.plot(lat_diff_list, std_list, "ko")
    plt.xlabel("Average latitudinal difference")
    plt.ylabel("Standard deviation of histograms")
    plt.title("stds over latitudinal distance")
    # plt.legend(["Linear regression", "data points"])

    plt.figure(2)
    plt.plot(long_diff_list, std_list, "ko")
    plt.xlabel("Average longitudinal difference")
    plt.ylabel("Standard deviation of histograms")
    plt.title("stds over longitudinal distance")

    plt.figure(3)
    plt.plot(dist_list, std_list, "ko")
    plt.xlabel("Great circle distance [m]")
    plt.ylabel("Standard deviation of histograms")
    plt.title("Stds over great circle distance")

    plt.figure(4)
    plt.plot(shift_list, std_list, "ko")
    plt.xlabel("Time between satellites [s]")
    plt.ylabel("std of histograms")
    plt.title("STD timeshift")

    plt.figure(5)
    plt.plot(shift_list, lat_diff_list, "ko")
    plt.title("Mean absolute latitudinal difference")
    plt.xlabel("Time difference $\Delta t$ [s]")
    plt.ylabel("$<|\Delta LAT|>$")
    plt.grid("on")
    plt.savefig("Figures/matfigs/latdist_over_timediff_all.pdf")
    
    plt.figure(6)
    plt.plot(shift_list, long_diff_list, "ko")
    plt.title("long diff")


    
    m, c, mdiv, cdiv = pro.linear_regression(shift_list, dist_list)
    xs = np.linspace(np.min(shift_list), np.max(shift_list), 100)
    print("Linear regression for distance over time, a = %g pm %g" % (m, mdiv))
    print("b = %g pm %g" % (c, cdiv))
    plt.figure(7)
    plt.plot(shift_list, dist_list, "ko")
    plt.plot(xs, xs*m + c)
    plt.title("Distance over time")
    plt.xlabel("Time difference $\Delta t$ [s]")
    plt.ylabel("Great circle distance [m]")
    plt.grid("on")
    plt.legend(["Average great circle distance", "Linear regression"])
    plt.savefig("Figures/matfigs/distance_over_timediff_all_a=%g_pm%g_b=%g_pm%g.pdf" % (m, mdiv, c, cdiv))
    plt.show()






    stop = time.time()
    print(stop-start)


sigma_plotter_mat()
