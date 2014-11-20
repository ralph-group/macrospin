"""
Helper functions for dealing with crystals
"""
import numpy as np

def miller_indices(lmn_vector):
    """ Returns an array of Miller indicies from the [lmn] unit vector """
    return (lmn_vector/np.abs(lmn_vector)).astype(int)

def miller_indices_display(lmn_array):
    """ Returns a string in LaTex format for Miller indices """
    result = '$['
    for i in lmn_array:
        if i < 0:
            result += r"\bar{%d}" % i
        else:
            result += "%d" % i
    return result + "]$"

def all_miller_indices(magnitude=3):
    """ Generates all Miller indices below a certain magnitude
    Based on: http://stackoverflow.com/a/25655090
    """
    span = np.linspace(-magnitude, magnitude, num=(2*magnitude+1))
    arrays = [span, span, span]
    shape = (len(x) for x in arrays)

    ix = np.indices(shape, dtype=int)
    ix = ix.reshape(len(arrays), -1).T

    for n, arr in enumerate(arrays):
        ix[:, n] = arrays[n][ix[:, n]]
    
    return ix

def angle_between_directions(origin, directions, rotation='left'):
    """ Returns the angles in degrees between the 
    """
    pass
