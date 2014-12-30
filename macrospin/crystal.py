"""
Helper functions for dealing with crystals
"""
import numpy as np
from __init__ import normalize

l = np.array([1, 0, 0])
m = np.array([0, 1, 0])
n = np.array([0, 0, 1])


def miller_indices(lmn_vector):
    """ Returns an array of Miller indicies from the [lmn] unit vector """
    return (lmn_vector/np.abs(lmn_vector)).astype(int)


def miller_indices_display(lmn_array):
    """ Returns a string in LaTex format for Miller indices """
    result = '$['
    for i in lmn_array:
        if i < 0:
            result += r"\bar{%d}" % abs(i)
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

    # Return only indices with sum greater than 0
    return ix[np.where(np.abs(ix).sum(axis=1) > 0)]


def angle_between_directions(start, normal, directions, rotation='right',
                             decimal_places=5):
    """ Returns the angles in degrees between the start the array of
    directions based on rotating a particular direction along a plane
    with a specified surface normal
    """
    orthogonal = np.cross(normal, start)
    orthogonal_product = np.round(directions.dot(orthogonal), decimal_places)
    angles = np.arccos(np.round(directions.dot(start),
                       decimal_places))*180/np.pi
    if rotation == 'right':
        subset = np.where(orthogonal_product > 0)
        angles[subset] = 360 - angles[subset]
    elif rotation == 'left':
        subset = np.where(orthogonal_product < 0)
        angles[subset] = 360 - angles[subset]
    else:
        raise ValueError('Inappropriate rotation direction specified')
    return angles


def in_plane_indices(normal, magnitude=1, decimal_places=5):
    """ Returns an array of Miller indices which lie in the
    plane specified by a surface normal vector. The magnitude
    specifies the maximum magnitude of the Miller indice searched.
    """
    indices = all_miller_indices(magnitude)
    return indices[np.where(
        np.round(normalize(indices).dot(normal), decimal_places) == 0
    )]


def in_plane_indices_and_angles(start, normal, rotation='right',
                                magnitude=1, decimal_places=5):
    """ Returns the Miller indices and angles in degrees for in-plane
    directions, based on a surface normal and a starting vector in the plane
    """
    indices = in_plane_indices(normal, magnitude, decimal_places)
    angles = angle_between_directions(start, normal, normalize(indices))
    mask = angles.argsort()
    return in_plane_indices[mask], angles[mask]
