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
    start_index = 1920 #index of 16 minutes
    stop_index = 3120 #index of 26 minutes

    classy = SWARMprocess()
    max_shift_ba = classy.timeshift(cdfB["Ne"][:100000], cdfA["Ne"][:100000], cdfA["Timestamp"][:100000],\
                        start = start_index, stop = stop_index, shifts = 5000 )
    max_shift_bc = classy.timeshift(cdfB["Ne"][:100000], cdfC["Ne"][:100000], cdfC["Timestamp"][:100000],\
                        start = start_index, stop = stop_index, shifts = 5000)

    plt.plot(cdfB["Timestamp"][start_index:stop_index], cdfB["Ne"][start_index:stop_index])
    plt.plot(cdfB["Timestamp"][start_index:stop_index],
            cdfA["Ne"][start_index + max_shift_ba: stop_index + max_shift_ba])
    plt.plot(cdfB["Timestamp"][start_index:stop_index],
            cdfC["Ne"][start_index + max_shift_bc:stop_index + max_shift_bc])
    plt.xlabel("time of satellite B")
    plt.ylabel("Electron Density")
    plt.title("Electron densities at north pole")
    plt.legend(["Satellite B", "Satellite A", "Satellite C"])
    plt.savefig("/home/aksel/Documents/Master/Spacemaster/Figures/polar_density.png")
    plt.show()

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
    classy = SWARMprocess()
    seconds = classy.stamp_to_sec(cdfA["Timestamp"][start:stop])

    lat_diff = latB - latA
    dist = classy.distance(latA, 0, radA, latB, 0, radA)

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
    plt.plot(seconds, lat_diff/np.max(lat_diff))
    plt.xlabel("seconds")
    plt.ylabel("Normalized distance and latitude")
    plt.legend(["Distance", "lat A", "lat B", "lat diff"])
    plt.show()
