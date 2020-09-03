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

M = int(50000)
N = int(5000)
NeA = np.array(cdfA["Ne"][:M])
NeB = np.array(cdfB["Ne"][:M])
times = cdfA["Timestamp"][:M]

classy = SWARMprocess()
shift_ba = classy.timeshift(NeB, NeA, times, start = 0, stop = 10000)


latA = np.array(cdfA["Latitude"][shift_ba:N + shift_ba])
longA = np.array(cdfA["Longitude"][shift_ba:N + shift_ba])
radA = np.array(cdfA["Radius"][shift_ba:N + shift_ba])

latB = np.array(cdfB["Latitude"][:N])
longB = np.array(cdfB["Longitude"][:N])
radB = np.array(cdfB["Radius"][:N])

NeA = NeA[shift_ba:N + shift_ba]
NeB = NeB[:N]

dist_ba = classy.distance(latB, longB, radB, latA, longA, radA)
seconds = classy.stamp_to_sec(times[:N])

# plt.plot(seconds, dist_ba/1e3)
# plt.xlabel("Time [s]")
# plt.ylabel("Distance [km]")
# plt.show()
#
# plt.plot(seconds, NeB)
# plt.plot(seconds, NeA)
# plt.xlabel("Time [s]")
# plt.ylabel("Ne")
# plt.show()
#
# plt.plot(seconds, radA)
# plt.plot(seconds, radB)
# plt.xlabel("Time [s]")
# plt.ylabel("Radius")
# plt.show()
#
# plt.plot(seconds, latA)
# plt.plot(seconds, latB)
# plt.xlabel("Time [s]")
# plt.ylabel("Latitude")
# plt.show()
#
# plt.plot(seconds, (longB - longA))
# plt.xlabel("Time [s]")
# plt.ylabel("relative Longitude")
# plt.show()
