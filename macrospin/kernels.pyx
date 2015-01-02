#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Kernel Cython classes for efficiently evolving the macrospin
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014-2015 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# distutils: language = c++
cimport cython
import numpy as np
cimport numpy as np

from macrospin.float3 cimport *
from macrospin cimport fields

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

    @property
    def t_sec(self):
        """ Returns the simulation time in seconds
        """
        return self.t/self.parameters['time_conversion']


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
                m = m + DT*(hxm + DAMPING*mxhxm)
            m.normalize()
            moments[i][0] = m.x
            moments[i][1] = m.y
            moments[i][2] = m.z

        self.m = np.array([m.x, m.y, m.z], dtype=np.float32)
        self.t += external_steps*internal_steps*DT

    def stabilize(self, const float threshold=1e-5, const long internal_steps=250):
        """ Runs the kernel until the cross product of the moment with
        the previous moment has a magnitude less than the threshold
        """

        cdef:
            long i, j
            long external_steps = 0
            float3 hxm, mxhxm
            float3 m = make_float3(self.m) # Initial orientation
            float3 h_ext = make_float3(self.parameters['Hext'])
            float DT = self.parameters['dt']
            float DAMPING = self.parameters['damping']

        m.normalize() # Initial normalization

        while True:
            mi = m
            for j in range(internal_steps):
                h_eff = h_ext
                hxm = h_eff.cross(m)
                mxhxm = m.cross(hxm)
                m = m + (hxm + mxhxm*DAMPING)*DT
            m.normalize()
            external_steps += 1
            if mi.cross(m).mag() <= threshold:
                break

        self.m = np.array([m.x, m.y, m.z], dtype=np.float32)
        self.t += external_steps*internal_steps*DT


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Anisotropy Kernel
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


class AnisotropyKernel(Kernel):


    def evolve(self, float[:,::1] moments, const long internal_steps=250):
        """ Fills an empty array of moment orientations with simulation results.
        Internal steps are taken for each data point, after which the moment is
        normalized.
        """
        
        cdef:
            long i, j
            long external_steps = moments.shape[0]
            float DT = self.parameters['dt']
            float DAMPING = self.parameters['damping']            
            float3 hxm, mxhxm, h_eff
            float3 m = make_float3(self.m) # Initial orientation
            float3 h_ext = make_float3(self.parameters['Hext'])
            float3 N = make_float3(self.parameters['Nd'])
            float hu1 = self.parameters['Hu1']
            float hu2 = self.parameters['Hu2']
            float3 eu = make_float3(self.parameters['eu'])
            float m_eu

        m.normalize() # Initial normalization

        for i in range(external_steps):
            for j in range(internal_steps):
                m_eu = m.dot(eu)
                h_eff = h_ext - m*N + eu*(m_eu*hu1) + eu*(m_eu*m_eu*m_eu*hu2)
                hxm = h_eff.cross(m)
                mxhxm = m.cross(hxm)
                m = m + (hxm + mxhxm*DAMPING)*DT
            m.normalize()
            moments[i][0] = m.x
            moments[i][1] = m.y
            moments[i][2] = m.z

        self.m = np.array([m.x, m.y, m.z], dtype=np.float32)
        self.t += external_steps*internal_steps*DT