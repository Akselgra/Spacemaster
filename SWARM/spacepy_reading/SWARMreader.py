"""
Currently only runs on Aksels laptop.
"""
import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
#setting path to cdf library
from getpass import getuser
usrname = getuser()
import os
os.environ["CDF_LIB"] = "/home/" + usrname +  "/Libraries/cdf/cdf36_3-dist/lib"
from spacepy import pycdf
from scipy.stats import pearsonr

data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"
# cdfA = pycdf.CDF("Data/Sat_A/SW_OPER_EFIA_LP_1B_20131202T101113_20131202T140109_0501.CDF/SW_OPER_EFIA_LP_1B_20131202T101113_20131202T140109_0501_MDR_EFI_LP.cdf")
# cdfB = pycdf.CDF("Data/Sswarm-diss.eo.esa.intat_B/SW_OPER_EFIB_LP_1B_20131202T114445_20131202T153609_0501.CDF/SW_OPER_EFIB_LP_1B_20131202T114445_20131202T153609_0501_MDR_EFI_LP.cdf")
# cdfC = pycdf.CDF("Data/Sat_C/SW_OPER_EFIC_LP_1B_20131204T094004_20131204T223759_0501.CDF/SW_OPER_EFIC_LP_1B_20131204T094004_20131204T223759_0501_MDR_EFI_LP.cdf")

cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIA_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIB_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"
cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501.CDF/SW_OPER_EFIC_LP_1B_20131221T000000_20131221T235959_0501_MDR_EFI_LP.cdf"

cdfA = pycdf.CDF(cdfA_path)
cdfB = pycdf.CDF(cdfB_path)
cdfC = pycdf.CDF(cdfC_path)

print(cdfA["U_orbit"].attrs)

def correlation_plotter(cdfA = pycdf.CDF(cdfA_path), cdfB = pycdf.CDF(cdfB_path), cdfC = pycdf.CDF(cdfC_path)):
    N = int(1e5)
    NeA = cdfA["Ne"][:N]
    NeB = cdfB["Ne"][:N]
    NeC = cdfC["Ne"][:N]
    time = cdfA["Timestamp"][:N]

    classy = SWARMprocess()
    corr_vec, shiftvec = classy.correlator(NeA, NeB, time)
    corr_vec_2, shiftvec_2 = classy.correlator(NeA, NeC, time)
    corr_vec_3, shiftvec_3 = classy.correlator(NeB, NeC, time)

    plt.plot(shiftvec/2, corr_vec)
    plt.plot(shiftvec/2, corr_vec_2)
    plt.plot(shiftvec/2, corr_vec_3)
    plt.xlabel("Timeshift [s]")
    plt.ylabel("Pearson R coefficient")
    plt.title("SWARM Ne correlation coefficients")
    plt.legend(["A and B", "A and C", "B and C"])
    plt.savefig("/home/aksel/Documents/Master/Spacemaster/Figures/correlations.png")
    plt.show()

def polar_plotter():
    """
    Shifts data to highest correlation and plots electron density at high latitudes.
    """
    start = 1920 #index of 16 minutes
    stop = 3120 #index of 26 minutes
    N = int(1e5)

    classy = SWARMprocess()
    ba_corrshift = classy.timeshift(cdfB["Ne"][:N], cdfA["Ne"][:N], cdfA["Timestamp"][:N],\
                        start = 0, stop = 5000, shifts = 5000 )
    bc_corrshift = classy.timeshift(cdfB["Ne"][:N], cdfC["Ne"][:N], cdfC["Timestamp"][:N],\
                        start = 0, stop = 5000, shifts = 5000)

    latA = cdfA["Latitude"][:N]
    latB = cdfB["Latitude"][:N]
    latC = cdfC["Latitude"][:N]

    ba_distshift = classy.timeshift_latitude(latB, latA, start, stop, shifts = 10000)
    bc_distshift = classy.timeshift_latitude(latB, latC, start, stop, shifts = 10000)

    NeA_corr = cdfA["Ne"][start+ba_corrshift:stop+ba_corrshift]
    NeC_corr = cdfC["Ne"][start+bc_corrshift:stop+bc_corrshift]

    NeA_dist = cdfA["Ne"][start+ba_distshift:stop+ba_distshift]
    NeC_dist = cdfC["Ne"][start+bc_distshift:stop+bc_distshift]

    NeB = cdfB["Ne"][start:stop]
    seconds = classy.stamp_to_sec(cdfA["Timestamp"][start:stop])

    plt.plot(latB[start:stop], NeB)
    plt.plot(latB[start:stop], NeA_corr)
    plt.plot(latB[start:stop], NeC_corr)
    plt.xlabel("Seconds since midnight 20.12.13")
    plt.ylabel("Electron density [cm⁻¹]")
    plt.legend(["Sat B", "Sat A", "Sat C"])
    plt.title("Sat A shifted %g indices forwards, sat C %g" %\
                                    (ba_corrshift, bc_corrshift))
    plt.show()

    plt.plot(seconds, NeB)
    plt.plot(seconds, NeA_dist)
    plt.plot(seconds, NeC_dist)
    plt.xlabel("Seconds since midnight 20.12.13")
    plt.ylabel("Electron density [cm⁻¹]")
    plt.title("Sat A shifted %g indices forwards, sat C %g" %\
                                  (ba_distshift, bc_distshift))
    plt.show()
    # plt.plot(cdfB["Timestamp"][start_index:stop_index], cdfB["Ne"][start_index:stop_index])
    # plt.plot(cdfB["Timestamp"][start_index:stop_index],
    #         cdfA["Ne"][start_index + max_shift_ba: stop_index + max_shift_ba])
    # plt.plot(cdfB["Timestamp"][start_index:stop_index],
    #         cdfC["Ne"][start_index + max_shift_bc:stop_index + max_shift_bc])
    # plt.xlabel("time of satellite B")
    # plt.ylabel("Electron Density")
    # plt.title("Electron densities at north pole")
    # plt.legend(["Satellite B", "Satellite A", "Satellite C"])
    # plt.savefig("/home/aksel/Documents/Master/Spacemaster/Figures/polar_density.png")
    # plt.show()

def distanceplotter():
    M = int(100000)
    N = int(50000)
    NeA = np.array(cdfA["Ne"][:M])
    NeB = np.array(cdfB["Ne"][:M])
    times = cdfA["Timestamp"][:M]

    classy = SWARMprocess()
    shift_ba = classy.timeshift(NeB, NeA, times, start = 0, stop = 10000)


    latA = np.array(cdfA["Latitude"][shift_ba:N+shift_ba])
    longA = np.array(cdfA["Longitude"][shift_ba:N+shift_ba])
    radA = np.array(cdfA["Radius"][shift_ba:N+shift_ba]) - classy.Re

    latB = np.array(cdfB["Latitude"][:N])
    longB = np.array(cdfB["Longitude"][:N])
    radB = np.array(cdfB["Radius"][:N]) - classy.Re

    NeA = NeA[shift_ba:N + shift_ba]
    NeB = NeB[:N]

    dist_ba = classy.distance(latB, longB, radB, latA, longA, radA)
    seconds = classy.stamp_to_sec(times[:N])

    plt.plot(seconds, dist_ba/1e3)
    plt.xlabel("Time [s]")
    plt.ylabel("Distance [km]")
    plt.title("Distance between closest measurement points for A and B")
    plt.savefig("/home/aksel/Documents/Master/Spacemaster/Figures/distance_time.png")
    plt.show()

    plt.plot(latB, dist_ba/1e3)
    plt.xlabel("Latitude")
    plt.ylabel("Distance [km]")
    plt.title("Distance between measurements for A and B over latitude")
    plt.savefig("/home/aksel/Documents/Master/Spacemaster/Figures/distance_latitude.png")
    plt.show()

def timediff_inspect():
    """
    Compares methods to find time difference between satellites
    """
    M = int(100000)
    NeA = np.array(cdfA["Ne"][:M])
    NeB = np.array(cdfB["Ne"][:M])
    NeC = np.array(cdfC["Ne"][:M])
    times = cdfA["Timestamp"][:M]

    classy = SWARMprocess()
    shift_ba = classy.timeshift(NeB, NeA, times)
    shift_bc = classy.timeshift(NeB, NeC, times)

    LatA = np.array(cdfA["Latitude"][:M])
    LatB = np.array(cdfB["Latitude"][:M])
    LatC = np.array(cdfC["Latitude"][:M])

    latshift_ba = classy.timeshift_latitude(LatB, LatA, shifts = 10000)
    latshift_bc = classy.timeshift_latitude(LatB, LatC, shifts = 10000)

    shifts = int(30000)
    meandist = np.zeros(shifts)
    indices = np.arange(shifts)
    for i in range(shifts):
        meandist[i] = np.mean(np.abs(LatB[0:40000] - LatA[i:40000 + i]))


    corrvec, shiftvec = classy.correlator(NeB, NeA, times)
    meandist = 1 - meandist/np.max(meandist)

    plt.plot(indices/2, meandist)
    plt.plot(shiftvec, corrvec)
    plt.xlabel("Seconds shifted")
    plt.ylabel("Normalized distance and correlation")
    plt.legend(["1 - normalized mean distance", "Pearson R coefficient"])
    plt.title("Correlation yields %g, distance yields %g" % (shift_ba, latshift_ba))
    plt.show()

def idek():
    """
    i dont even know
    """

    N = int(1e4)
    start = 1920
    stop = 3120
    shift = 90
    latA = cdfA["Latitude"][start+shift:stop+shift]
    latB = cdfB["Latitude"][start:stop]
    radA = cdfA["Radius"][start:stop]
    longA = cdfA["Longitude"][start+shift:stop+shift]
    longB = cdfB["Longitude"][start:stop]
    classy = SWARMprocess()
    seconds = classy.stamp_to_sec(cdfA["Timestamp"][start:stop])

    lat_diff = latB - latA
    dist = classy.distance(latA, longA, radA, latB, longB, radA)

    latA1 = latA[:-1]
    latA2 = latA[1:]
    dt = seconds[1] - seconds[0]
    latAdiff = latA2 - latA1
    deri_latA = latAdiff/dt

    latB1 = latB[:-1]
    latB2 = latB[1:]
    dt = seconds[1] - seconds[0]
    latBdiff = latB2 - latB1
    deri_latB = latBdiff/dt

    secondish = seconds[:-1]

    deri_diff = deri_latB - deri_latA
    longdiff = longB - longA

    plt.plot(secondish, deri_latA)
    plt.plot(secondish, deri_latB)
    plt.plot(secondish, deri_diff)
    plt.ylabel("Time derivative of latitude")
    plt.xlabel("Seconds")
    plt.legend(["sat A", "sat B", "difference"])
    plt.show()

    plt.plot(seconds, dist/np.max(dist))
    plt.plot(seconds, latA/np.max(latA))
    plt.plot(seconds, latB/np.max(latB))
    plt.plot(seconds, lat_diff/np.max(lat_diff), ".")
    plt.plot(seconds, longdiff/np.max(longdiff))
    plt.xlabel("seconds")
    plt.ylabel("Normalized distance and latitude")
    plt.legend(["Distance", "lat A", "lat B", "lat diff", "long diff"])
    plt.show()

def funky():
    """
    Integrated fourier at obvious spatial difference
    """
    start = 1100*2 #index of 16 minutes
    stop = 1230*2 #index of 26 minutes
    N = int(1e5)

    classy = SWARMprocess()
    latA = cdfA["Latitude"][:N]
    latB = cdfB["Latitude"][:N]
    latC = cdfC["Latitude"][:N]

    ba_distshift = classy.timeshift_latitude(latB, latA, start, stop, shifts = 10000)
    bc_distshift = classy.timeshift_latitude(latB, latC, start, stop, shifts = 10000)

    NeA = cdfA["Ne"][start+ba_distshift:stop+ba_distshift]
    NeC = cdfC["Ne"][start+bc_distshift:stop+bc_distshift]

    NeB = cdfB["Ne"][start:stop]
    seconds = classy.stamp_to_sec(cdfA["Timestamp"][start:stop])
    plt.plot(seconds, NeB/np.max(NeB), "b")
    plt.plot(seconds, NeA/np.max(NeA), "g")
    plt.plot(seconds, NeC/np.max(NeC), "r")

    n = 85
    minfreq = 0
    maxfreq = 1/3
    times, timyA = classy.fft_time_integral(NeA, n, 2, minfreq = minfreq, maxfreq = maxfreq)
    times, timyB = classy.fft_time_integral(NeB, n, 2, minfreq = minfreq, maxfreq = maxfreq)
    times, timyC = classy.fft_time_integral(NeC, n, 2, minfreq = minfreq, maxfreq = maxfreq)
    times += start/2

    diffBA = classy.relative_diff(timyB, timyA, abs = False)
    diffBC = classy.relative_diff(timyB, timyC, abs = False)
    diffAC = classy.relative_diff(timyA, timyC, abs = False)

    plt.plot(times, diffBA, "b.")
    plt.plot(times, diffBC, "g.")
    plt.plot(times, diffAC, "r.")
    plt.show()


def funky_fftint():
    """
    Specific integrated fft
    """
    start = 1165*2 #index of 16 minutes
    stop = 1185*2 #index of 26 minutes
    N = int(1e5)

    classy = SWARMprocess()
    latA = cdfA["Latitude"][:N]
    latB = cdfB["Latitude"][:N]
    latC = cdfC["Latitude"][:N]

    ba_distshift = classy.timeshift_latitude(latB, latA, start, stop, shifts = 10000)
    bc_distshift = classy.timeshift_latitude(latB, latC, start, stop, shifts = 10000)

    NeA = cdfA["Ne"][start+ba_distshift:stop+ba_distshift]
    NeC = cdfC["Ne"][start+bc_distshift:stop+bc_distshift]

    NeB = cdfB["Ne"][start:stop]
    seconds = classy.stamp_to_sec(cdfA["Timestamp"][start:stop])

    fs = 2
    n = len(NeA)

    fftA = np.fft.fft(NeA)[:int(n/2)]
    fftB = np.fft.fft(NeB)[:int(n/2)]
    fftC = np.fft.fft(NeC)[:int(n/2)]
    omegas = np.linspace(-1, 1, n)[:int(n/2)]
    fftA = np.abs(fftA)
    fftB = np.abs(fftB)
    fftC = np.abs(fftC)

    df = 0.1
    deltaf = 0.1

    nf = int((1- deltaf)/df)
    f0s = np.linspace(0, 0.9, nf)
    f1s = np.linspace(df, 1, nf)

    A_ints = np.zeros_like(f0s)
    B_ints = np.zeros_like(f0s)
    C_ints = np.zeros_like(f0s)

    d_omega = omegas[1] - omegas[0]
    for i in range(len(f0s)):
        minfreq = f0s[i]
        maxfreq = f1s[i]
        n_minfreq = int(minfreq/d_omega)
        n_maxfreq = int(maxfreq/d_omega)
        A_ints[i] = np.sum(fftA[n_minfreq:n_maxfreq])*d_omega
        B_ints[i] = np.sum(fftB[n_minfreq:n_maxfreq])*d_omega
        C_ints[i] = np.sum(fftC[n_minfreq:n_maxfreq])*d_omega


    plt.figure(0)
    plt.plot(seconds, NeB)
    plt.plot(seconds, NeA)
    plt.plot(seconds, NeC)
    plt.xlabel("Time of sat B [s]")
    plt.ylabel("Electron density [cm⁻¹]")
    plt.title("SWARM EFI LP data")

    plt.figure(1)
    plt.plot(f0s + deltaf/2, A_ints, "-o")
    plt.plot(f0s + deltaf/2, B_ints, "-o")
    plt.plot(f0s + deltaf/2, C_ints, "-o")
    plt.legend(["A", "B", "C"])
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Integrated fourier")
    plt.title("Integrated fourier at spatial difference")


    BAdiff = classy.relative_diff(B_ints, A_ints, abs = False, norm = 1)
    BCdiff = classy.relative_diff(B_ints, C_ints, abs = False, norm = 1)
    ACdiff = classy.relative_diff(A_ints, C_ints, abs = False, norm = 1)

    plt.figure(2)
    plt.axis([0, 1, -1, 1])
    plt.plot(f0s + deltaf/2, BAdiff, "-o")
    plt.plot(f0s + deltaf/2, BCdiff, "-o")
    plt.plot(f0s + deltaf/2, ACdiff, "-o")
    plt.legend(["BAdiff", "BCdiff", "ACdiff"])
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("normalized difference in integrated fourier")
    plt.title("Difference in integrated fourier at spatial difference")
    plt.show()



#funky_fftint()
