from SWARMprocess import SWARMprocess
import numpy as np


def test_distance():
    """
    test function for distance method in SAWRMprocess.
    """
    r = 2
    lat = 90
    long1 = 90
    long2 = 270

    analdist = 4

    classy = SWARMprocess()
    compdist = classy.distance(lat, long1, r, lat, long2, r)

    epsilon = 1e-4
    testy = (np.abs(analdist - compdist) < epsilon)
    msg = "analytical distance not equal to calculated distance.\
     Calculated %g, expected %g" % (compdist, analdist)
    assert testy, msg

def test_maxdiff():
    """
    test function for maxdiff method in SAWRMprocess
    """

    array1 = np.array([1,2,3,2,1])
    array2 = np.array([1,2,3,4,5])
    diff = -2
    object = SWARMprocess()
    calcdiff = object.maxdiff(array1, array2)[1]
    eps = 1e-5
    testy = np.abs(calcdiff - diff) < eps
    msg = "Calculated difference is %g while true dist is %g" % (calcdiff, diff)
    assert testy, msg

print("Running tests")
test_distance()
test_maxdiff()
print("Tests passed")
