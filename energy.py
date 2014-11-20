"""
Functions for the free energy of a macrospin magnet
"""
from __future__ import division
import numpy as np

def uniaxial_anisotropy(m, u=np.array([0,0,1]), Ku1=0, Ku2=0):
    """ Returns the uniaxial anisotropy part of the free energy in erg/cc
    m: Reduced magnetization M/Ms (dimensionless)
    u: Uniaxial direction (unit vector)
    Ku1: Uniaxial anisotropy energy 1 (erg/cc)
    Ku2: Uniaxial anisotropy energy 2 (erg/cc)    
    """
    Fu = -Ku1*(u.dot(m)**2)
    Fu += -Ku2*(u.dot(m)**4)
    return Fu
    
def cubic_anisotropy(m, c1, c2, Kc1=0, Kc2=0, Kc3=0):
    """ Returns the cubic anisotropy part of the free energy in erg/cc
    m: Reduced magnetization M/Ms (dimensionless)
    c1: Orthogonal cubic direction 1 (unit vector)
    c2: Orthogonal cubic direction 2 (unit vector)
    c3: Orthogonal cubic direction 3 (unit vector)
    Kc1: Cubic anisotropy energy 1 (erg/cc)
    Kc2: Cubic anisotropy energy 2 (erg/cc) 
    Kc2: Cubic anisotropy energy 3 (erg/cc)   
    """
    c3 = np.cross(c1, c2)
    
    mc1 = c1.dot(m)
    mc2 = c2.dot(m)
    mc3 = c3.dot(m)
    
    Fc =  Kc1*((mc1**2)*(mc2**2))
    Fc += Kc1*((mc2**2)*(mc3**2))
    Fc += Kc1*((mc1**2)*(mc3**2))
    
    Fc += Kc2*((mc1**2)*(mc2**2)*(mc3**2))
    
    Fc =  Kc3*((mc1**4)*(mc2**4))
    Fc += Kc3*((mc2**4)*(mc3**4))
    Fc += Kc3*((mc1**4)*(mc3**4))
    
    return Fc
