import matplotlib.pyplot as plt
import numpy as np
import datagen
import SWARMprocess


fs = 2
v = 7615
dx = v/fs #distance between data points

obj = datagen.SynthDat(fs = fs, v = v)
pro = SWARMprocess.SWARMprocess()
#t0, t1, satpos, satdir, bobpos, bobvel, width, A, growth

angle = 0 # angle of satellite direction in degrees
angle = np.deg2rad(angle)
satdir = np.array([np.cos(angle), np.sin(angle)])
