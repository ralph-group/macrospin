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

from macrospin.types cimport * # includes float3 and its operators
from macrospin cimport field, torque, energy, solvers


cdef struct Step:
    """ Stores the results of a time step operation

    m - moment unit vector
    t - simulation time
    torque - torque at (t, m)
    eps - error

    """
    float3 m, torque
    float t, eps


cdef class Kernel(object):
    """ Encapsulates the time evolution algorithm for solving the
    Landau-Liftshitz equation
    """
    cdef:
        object parameters
        Step previous, current
        void *step_func(Kernel*)
        float dt


    def __cinit__(self, object parameters, char* step_method='RK23'):
        self.parameters = parameters.normalize()
        self._load()

        self.set_step_func(step_method)

        self.reset()


    # TODO: def __dealloc__(self):


    def _load(self):
        """ Loads the parameters from the kernel (to be extended)
        """
        self.dt = self.parameters['dt']


    def run(self, time=None, internal_steps=250):
        """ Run the simulation for a given time
        """
        time *= self.parameters['time_conversion']
        steps = int(time/self.parameters['dt']/internal_steps)

        cdef float[:,::1] moments = np.zeros((steps, 3), dtype=np.float32)
        times = self.times(moments, internal_steps)

        self._evolve(&moments[0][0], steps, internal_steps)

        return times, moments


    cdef void _evolve(self, float* moments_ptr, long steps, long internal_steps):
        """ Takes steps

        **moments - pointer to array of array of floats
        external_steps -
        internal_steps -

        """

        for i in range(external_steps):
            for j in range(internal_steps):
                self.previous = self.current
                self.step_func(self)

            moments_ptr[3*i] = self.current.m.x;
            moments_ptr[3*i+1] = self.current.m.y;
            moments_ptr[3*i+2] = self.current.m.z;             


    cdef void set_step_func(step_func_name):
        """ Sets the step function based on its name """
        if step_func_name == "Euler":
            self.step_func = solvers.euler_step
        elif step_func_name == "Huen":
            self.step_func = solvers.huen_step
        elif step_func_name == "RK23":
            self.step_func = solvers.rk23_step
        elif step_func_name == "RK4":
            self.step_func = solver.rk4_step


    def times(self, moments, internal_steps=250):
        """ Returns an array of times that correspond to the moments array
        that is about to be run
        """
        return (self.t + internal_steps*self.parameters['dt']*(
                1 + np.arange(moments.shape[0])))/(
                self.parameters['time_conversion'])


    cdef float3 field(float t, float3 m):
        float3 field
        return field


    cdef float3 torque(float t, float3 m):
        float3 field
        return field    


    cdef float energy(float t, float3 m):
        return energy.zeeman(m, self.field(t, m))


    @property
    def m(self):
        """ Returns the moment unit vector as a numpy array
        """
        cdef float3 m = self.current.m
        return np.array([m.x, m.y, m.z], dtype=np.float32)


    @property
    def t(self):
        """ Returns the simulation time in units of (gamma Ms)
        """
        return self.current.t


    @property
    def t_sec(self):
        """ Returns the simulation time in seconds
        """
        return self.current.t/self.parameters['time_conversion']


    def reset(self):
        """ Resets the kernel to the initial conditions
        """
        self.current.m = make_float3(self.parameters['m0'])
        self.current.t = 0.0
        self.current.torque = self.torque(self.current.t, self.current.m)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Basic Kernel
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


cdef class BasicKernel(Kernel):
    cdef:
        float alpha
        float3 hext, Nd


    def _load(self):
        """ Loads the parameters from the kernel (overwritting Kernel.load)
        """
        self.dt = self.parameters['dt']
        self.alpha = self.parameters['damping']
        self.hext = make_float3(self.parameters['Hext'])
        self.Nd = make_float3(self.parameters['Nd'])


    cdef float3 field(float t, float3 m):
        return self.hext + field.demagnetization(m, self.Nd)


    cdef float3 torque(float t, float3 m):
        cdef float3 heff = self.field(t, m)
        return torque.landau_lifshitz(m, heff, self.alpha)