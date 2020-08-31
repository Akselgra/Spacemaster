"""
Currently only runs on Aksels laptop.
"""
import numpy as np
import matplotlib.pyplot as plt
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

print(cdfA)
