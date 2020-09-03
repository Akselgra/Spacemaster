from SWARMprocess import SWARMprocess
import numpy as np


def test_distance():
    """
    test function for distance method in SAWRMprocess.
    """
    r = 2
    lat = 90
    long1 = 0
    long2 = 180

    analdist = 4

    classy = SWARMprocess()
    compdist = classy.distance(lat, long1, r, lat, long2, r)

    epsilon = 1e-4
    testy = (np.abs(analdist - compdist) < epsilon)
    msg = "analytical distance not equal to calculated distance"
    assert testy, msg

if __name__ == "__main__":
    test_distance()
