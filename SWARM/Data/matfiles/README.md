Data files from SWARM.\n
.mat files containing data from the ESA SWARM mission.\n
Each .mat file contains one day of data for all 3 SWARM satellites.\n
Data is obtained by processing original CDF files. Geomagnetic coordinates are obtained using the aacgmv2 library.\n
Data includes:\n
NeA, NeB, NeC; - Electron density data from EFI\n
longA, longB, longC; - geocentric longitude coordinates in ITRF\n
latA, latB, latC; - geocentric latitude coordinates in ITRF\n
radA, radB, radC; - distance from earth's center in ITRF\n
velA, velB, velC; - orbital velocity of satellites\n
altA, altB, altC; - altitude of satellites (distance from sea level)\n
mlatA, mlatB, mlatC; - geomagnetic latitude coordinates\n
mlongA, mlongB, mlongC; - geomagnetic longitude coordinates\n
mltA, mltB, mltC; - magnetic local time\n
secondsA, secondsB, secondsC; - Time of measurement in seconds since midnight UTC.\n
