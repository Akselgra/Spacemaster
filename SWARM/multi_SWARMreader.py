"""
CDF_LIB must be in the correct file path to run
"""
import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess
from general_SWARMreader import GenSWARMread
#setting path to cdf library
from getpass import getuser
usrname = getuser()
import os
os.environ["CDF_LIB"] = "/home/" + usrname +  "/Libraries/cdf/cdf36_3-dist/lib"
from spacepy import pycdf


class MultiSWARM():
    """
    Class for reading multiple files of SWARM data.
    Designed for specific file structures.
    """

    def __init__(self, year1, month1, day1, year2, month2, day2):
        """
        parameters:
            year1 - int; year of first data file
            month1 - int; month of first data file
            day1 - int; day of first data file
            year2 - int; year of last data file
            month2 - int; month of last data file
            day2 - int; day of last data file
        """
        self.year1 = year1
        self.month1 = month1
        self.day1 = day1
        self.year2 = year2
        self.month2 = month2
        self.day2 = day2

        self.data_path = "/home/" + usrname +  "/Documents/Master/Swarm_Data"



    def filenames(self, year, month, day):
        """
        Takes a date and returns a list of file path strings.
        parameters:
            year - int; year of data files yyyy
            month - int; month of data files mm
            day - int; day of data files dd
        returns:
            files - list; file paths on the form ["satA", "satB", "satC"]
        """
        data_path = self.data_path

        year = str(year)
        if len(year) != 4:
            raise ValueError("year must be on form yyyy")

        month = str(month)
        if len(month) == 1:
            month = "0" + month

        day = str(day)
        if len(day) == 1:
            day = "0" + day

        date = year + month + day
        t0 = "T000000_"
        t1 = "T235959_0501"

        cdfA_path = data_path + "/Sat_A/SW_OPER_EFIA_LP_1B_" + date + t0 + date\
        + t1 + ".CDF/SW_OPER_EFIA_LP_1B_" + date + t0 + date + t1 + "_MDR_EFI_LP.cdf"

        cdfB_path = data_path + "/Sat_B/SW_OPER_EFIB_LP_1B_" + date + t0 + date\
         + t1 + ".CDF/SW_OPER_EFIB_LP_1B_" + date + t0 + date + t1 + "_MDR_EFI_LP.cdf"

        cdfC_path = data_path + "/Sat_C/SW_OPER_EFIC_LP_1B_" + date + t0 + date\
         + t1 + ".CDF/SW_OPER_EFIC_LP_1B_" + date + t0 + date + t1 + "_MDR_EFI_LP.cdf"

        files = [cdfA_path, cdfB_path, cdfC_path]

        return(files)

    def gen_filenames(self, ind):
        """
        Generalized version of filenames. Takes an index n and returns
        filepaths of sats A, B and C for the day of the given index in directory.
        example: ind 0 gives first day. ind 21 gives 21st day.
        parameters:
            ind - int; index of day in directory
        returns:
            files - list; list of filepaths on the form [satA, satB, satC]
        """

        pathA = self.data_path + "/Sat_A"
        dirsA = os.listdir(pathA)
        dirsA.sort()
        pathA = pathA + "/" + dirsA[ind]
        cdfA = os.listdir(pathA)
        cdfA.sort()
        pathA = pathA + "/" + cdfA[1]

        pathB = self.data_path + "/Sat_B"
        dirsB = os.listdir(pathB)
        dirsB.sort()
        pathB = pathB + "/" + dirsB[ind]
        cdfB = os.listdir(pathB)
        cdfB.sort()
        pathB = pathB + "/" + cdfB[1]

        pathC = self.data_path + "/Sat_A"
        dirsC = os.listdir(pathC)
        dirsC.sort()
        pathC = pathC + "/" + dirsC[ind]
        cdfC = os.listdir(pathC)
        cdfC.sort()
        pathC = pathC + "/" + cdfC[1]

        files = [pathA, pathB, pathC]

        print(pathA)
        print(pathB)
        print(pathC)
        return(files)

if __name__ == "__main__":
    object = MultiSWARM(1,1,1,2,2,2)
    files = object.gen_filenames(25)

    genread = GenSWARMread(files[0], files[1], files[2])
    hists, bins = genread.histmake()
    plt.bar(bins[0], hists[0], width = (bins[0][1] - bins[0][0])*0.9)
    plt.show()
