Data files from SWARM.\n
.mat files containing data from the ESA SWARM mission.\
Each .mat file contains one day of data for all 3 SWARM satellites.\
Data is obtained by processing original CDF files. Geomagnetic coordinates are obtained using the aacgmv2 library.\
Data can be retried using scipy.io.loadmat(filename)["attribute"][0]\
Data includes:
* NeA, NeB, NeC; - Electron density data from EFI
* longA, longB, longC; - geocentric longitude coordinates in ITRF
* latA, latB, latC; - geocentric latitude coordinates in ITRF
* radA, radB, radC; - distance from earth's center in ITRF
* velA, velB, velC; - orbital velocity of satellites
* altA, altB, altC; - altitude of satellites (distance from sea level)
* mlatA, mlatB, mlatC; - geomagnetic latitude coordinates
* mlongA, mlongB, mlongC; - geomagnetic longitude coordinates
* mltA, mltB, mltC; - magnetic local time
* secondsA, secondsB, secondsC; - Time of measurement in seconds since midnight UTC.
* BA_shift; - Time difference between satellites B and A
* BC_shift; - Time difference between satellites B and C
