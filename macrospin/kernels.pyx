#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Kernel Cython classes for efficiently evolving the macrospin
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# distutils: language = c++
cimport cython
import numpy as np
cimport numpy as np
from libcpp cimport queue


cdef extern from "float3.h":
    cdef cppclass float3:
        float3(x, y, z)
        float3()
        float x, y, z
        float3 operator+(float3)
        float3 operator-(float3)
        float3 operator*(float3)
        float3 operator*(float)
        float3 operator/(float)
        float3 cross(float3)
        float dot(float3)
        float mag()
        void normalize()


@cython.boundscheck(False)
@cython.wraparound(False)
cdef float3 make_float3(float[::1] l):
    cdef float3 r
    r.x, r.y, r.z = l[0], l[1], l[2]
    return r


class Kernel(object):
    """ Encapsulates the time evolution algorithm for solving the
    Landau-Liftshitz equation
    """

    def __init__(self, parameters):
        self.parameters = parameters.normalize()
        self.reset()


    def run(self, time=None, internal_steps=250):
        """ Run the simulation for a given time
        """
        time *= self.parameters['time_conversion']
        steps = int(time/self.parameters['dt']/internal_steps)

        moments = np.zeros((steps, 3), dtype=np.float32)
        times = self.times(moments, internal_steps)

        self.evolve(moments, internal_steps)

        return times, moments

    def times(self, moments, internal_steps=250):
        """ Returns an array of times that correspond to the moments array
        that is about to be run
        """
        return (self.t + internal_steps*self.parameters['dt']*(
                1 + np.arange(moments.shape[0])))/(
                self.parameters['time_conversion'])


    def reset(self):
        """ Resets the kernel to the initial conditions
        """
        self.m = self.parameters['m0']
        self.t = 0.0


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Basic Kernel
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


class BasicKernel(Kernel):


    def evolve(self, float[:,::1] moments, const long internal_steps=250):
        """ Fills an empty array of moment orientations with simulation results.
        Internal steps are taken for each data point, after which the moment is
        normalized.
        """
        
        cdef:
            long i, j
            long external_steps = moments.shape[0]
            float3 hxm, mxhxm
            float3 m = make_float3(self.m) # Initial orientation
            float3 h_ext = make_float3(self.parameters['Hext'])
            float DT = self.parameters['dt']
            float DAMPING = self.parameters['damping']

        m.normalize() # Initial normalization

        for i in range(external_steps):
            for j in range(internal_steps):
                h_eff = h_ext
                hxm = h_eff.cross(m)
                mxhxm = m.cross(hxm)
                m = m + (hxm + mxhxm*DAMPING)*DT
            m.normalize()
            moments[i][0] = m.x
            moments[i][1] = m.y
            moments[i][2] = m.z

        self.m = np.array([m.x, m.y, m.z], dtype=np.float32)
        self.t += external_steps*internal_steps*DT


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Anisotropy Kernel
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>