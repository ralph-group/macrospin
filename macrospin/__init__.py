"""
Includes helper functions for macrospin
"""
from __future__ import division
import os
import numpy as np
from mpl_toolkits.mplot3d import proj3d

__version__ = 0.1


def get_include():
    import macrospin
    return os.path.dirname(macrospin.__file__)+"/" # TODO: Eliminate +"/" hack
    

def brillouin(x, J):
    """ Returns the Brillouin function evaluated at x and J """
    f1 = (2*J+1)/(2*J)
    f2 = 1/(2*J)
    return f1/np.tanh(f1*x)-f2/np.tanh(f2*x)

def langevin(x):
    """ Returns the Langevin function evaluated at x and J """
    return np.coth(x)-1/x


def rotation_array(R_function, angles):
    """ Returns an array of rotation matrixes in correct order based on
    a rotation function and an array of angles in degrees

    """
    rotations = np.empty((len(angles), 3, 3))
    for i, angle in enumerate(angles):
        rotations[i] = R_function(angle)
    return rotations


def Rx(angle_degrees):
    """ Constructs a rotation matrix around x based on an angle in degrees """
    angle = angle_degrees/180*np.pi
    return np.matrix([
        [1, 0, 0],
        [0, np.cos(angle), -np.sin(angle)],
        [0, np.sin(angle), np.cos(angle)],
    ])


def Ry(angle_degrees):
    """ Constructs a rotation matrix around y based on an angle in degrees """
    angle = angle_degrees/180*np.pi
    return np.matrix([
        [np.cos(angle), 0, np.sin(angle)],
        [0, 1, 0],
        [-np.sin(angle), 0, np.cos(angle)],
    ])


def Rz(angle_degrees):
    """ Constructs a rotation matrix around z based on an angle in degrees """
    angle = angle_degrees/180*np.pi
    return np.matrix([
        [np.cos(angle), -np.sin(angle), 0],
        [np.sin(angle), np.cos(angle), 0],
        [0, 0, 1],
    ])


def normalize(array):
    """ Normalizes an array row-wise """
    if len(array.shape) == 1:
        return array/np.linalg.norm(array).T
    else:
        normals = np.linalg.norm(array, axis=1)
        normals = np.array(normals.reshape(normals.size, 1))
        return array / normals

