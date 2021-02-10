from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
from SWARMprocess import SWARMprocess

class MatReader(SWARMprocess):
    """
    class to read .mat files containing SWARM data.
    """

    def __init__(self, filename):
        """
        Parameters:
            filename:
                string; name of .mat file to be read
        """

        infile = loadmat(filename)
        self.NeA = infile["NeA"][0]
        self.NeB = infile["NeB"][0]
        self.NeC = infile["NeC"][0]

        self.longA = infile["longA"][0]
        self.longB = infile["longB"][0]
        self.longC = infile["longC"][0]

        self.latA = infile["latA"][0]
        self.latB = infile["latB"][0]
        self.latC = infile["latC"][0]

        self.radA = infile["radA"][0]
        self.radB = infile["radB"][0]
        self.radC = infile["radC"][0]

        self.velA = infile["velA"][0]
        self.velB = infile["velB"][0]
        #self.velC = infile["velC"][0]

        self.altA = infile["altA"][0]
        self.altB = infile["altB"][0]
        self.altC = infile["altC"][0]

        self.mlatA = infile["mlatA"][0]
        self.mlatB = infile["mlatB"][0]
        self.mlatC = infile["mlatC"][0]

        self.mlongA = infile["mlongA"][0]
        self.mlongB = infile["mlongB"][0]
        self.mlongC = infile["mlongC"][0]

        self.mltA = infile["mltA"][0]
        self.mltB = infile["mltB"][0]
        self.mltC = infile["mltC"][0]

        self.secondsA = infile["secondsA"][0]
        self.secondsB = infile["secondsB"][0]
        self.secondsC = infile["secondsC"][0]

        self.fs = 2

        self.BA_shift = self.timeshift_latitude(self.latB, self.latA)
        self.BC_shift = self.timeshift_latitude(self.latB, self.latC)

        print(len(self.latA))
        print(len(self.latC))

        print(self.BA_shift)
        print(self.BC_shift)

if __name__ == "__main__":
    file = "Data/matfiles/20131220.mat"
    object = MatReader(file)
