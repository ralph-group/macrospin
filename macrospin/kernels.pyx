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

from macrospin.cpp.c_types cimport * # includes float3 and its operators
from macrospin.cpp cimport c_kernels
from macrospin.cpp cimport c_solvers


cdef class Kernel(object):
    """ Encapsulates the time evolution algorithm for solving the
    Landau-Liftshitz equation
    """
    cdef:
        c_kernels.Kernel *kernel
        object parameters

    def __cinit__(self, object parameters, char* step_method='RK23'):
        self.parameters = parameters.normalize()
        self._load()

        self.kernel.set_step_func(step_method)

        self.reset()

    # TODO: def __dealloc__(self):

    def _load(self):
        """ Loads the parameters from the kernel (to be overwritten)
        """
        self.kernel = new c_kernels.Kernel()
        self.kernel.dt = self.parameters['dt']

    def run(self, time=None, internal_steps=250):
        """ Run the simulation for a given time
        """
        time *= self.parameters['time_conversion']
        steps = int(time/self.parameters['dt']/internal_steps)

        cdef float[:,::1] moments = np.zeros((steps, 3), dtype=np.float32)
        times = self.times(moments, internal_steps)

        self._evolve(&moments[0][0], steps, internal_steps)

        return times, moments

    cdef _evolve(self, float* moments_ptr, long steps, long internal_steps):
        self.kernel.evolve(moments_ptr, steps, internal_steps)

    def times(self, moments, internal_steps=250):
        """ Returns an array of times that correspond to the moments array
        that is about to be run
        """
        return (self.t + internal_steps*self.parameters['dt']*(
                1 + np.arange(moments.shape[0])))/(
                self.parameters['time_conversion'])

    @property
    def m(self):
        """ Returns the moment unit vector as a numpy array
        """
        cdef float3 m = self.kernel.current.m
        return np.array([m.x, m.y, m.z], dtype=np.float32)

    @property
    def t(self):
        """ Returns the simulation time in units of (gamma Ms)
        """
        return self.kernel.current.t

    @property
    def t_sec(self):
        """ Returns the simulation time in seconds
        """
        return self.kernel.current.t/self.parameters['time_conversion']


    def reset(self):
        """ Resets the kernel to the initial conditions
        """
        self.kernel.current.m = make_float3(self.parameters['m0'])
        self.kernel.current.t = 0.0
        self.kernel.current.torque = self.kernel.torque(self.kernel.current.t, 
            self.kernel.current.m)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Basic Kernel
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


cdef class BasicKernel(Kernel):

    def _load(self):
        """ Loads the parameters from the kernel (overwritting Kernel.load)
        """
        self.kernel = <c_kernels.Kernel *> new c_kernels.BasicKernel()
        self.kernel.dt = self.parameters['dt']
        (<c_kernels.BasicKernel *>self.kernel).alpha = self.parameters['damping']
        (<c_kernels.BasicKernel *>self.kernel).hext = make_float3(self.parameters['Hext'])
        (<c_kernels.BasicKernel *>self.kernel).Nd = make_float3(self.parameters['Nd'])

    cdef _evolve(self, float* moments_ptr, long steps, long internal_steps):
        (<c_kernels.BasicKernel *>self.kernel).evolve(moments_ptr, steps, internal_steps)