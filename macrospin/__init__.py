"""
Includes helper functions for macrospin
"""
import numpy as np

def rotation_array(R_function, angles):
    """ Returns an array of rotation matrixes in correct order based on
    a rotation function and an array of angles in degrees
    """
    rotations = np.empty((len(angles),3,3))
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
