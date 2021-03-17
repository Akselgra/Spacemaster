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
    
    
def test_great_circle_distance():
    """
    test function for great circle distance method in SWARMprocess
    """
    object = SWARMprocess()
    r = object.Re #earth radius
    lat1 = -2.5856855; long1 = -18.1962077
    lat2 = -6.405082999999999; long2 = -18.0331106
    
    computed = object.great_circle_distance(lat1, long1, r, lat2, long2, r)/1000
    expected = 425.08214222811756
    
    eps = 1
    testy = np.abs(computed - expected) < eps
    msg = "calculated great circle distance is not what we expected"
    
    print(computed)
    print(expected)
    assert testy, msg
    

print("Running tests")
test_distance()
test_maxdiff()
test_great_circle_distance()
print("Tests passed")
