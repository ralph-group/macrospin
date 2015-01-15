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
from macrospin.parameters import NormalizedParameters

from libc.math cimport sin, cos


cdef class Kernel:
    """ Encapsulates the time evolution algorithm for solving the
    Landau-Liftshitz equation
    """

    def __init__(self, object parameters, step_method="RK23"):
        self.raw_parameters = parameters
        self.parameters = NormalizedParameters(parameters)
        self._load()

        self.set_step_func(step_method)

        self.reset()


    # TODO: def __dealloc__(self):

    def set_step_func(self, name):
        if name == "Euler": self.step_func = &(solvers.euler_step)
        elif name == "Huen": self.step_func = &(solvers.huen_step)
        elif name == "RK23": self.step_func = &(solvers.rk23_step)
        elif name == "RK4": self.step_func = &(solvers.rk4_step)
        else: self.step_func = &(solvers.rk23_step) # Default to RK23


    def _load(self):
        """ Loads the parameters from the kernel (to be extended)
        """
        self.dt = self.parameters['dt']


    def reset(self):
        """ Resets the kernel to the initial conditions
        """
        self.current.m = make_float3(self.parameters['m0'])
        self.current.m.normalize()
        self.current.t = 0.0
        self.current.torque = self.torque(self.current.t, self.current.m)


    def run(self, time=None, internal_steps=250):
        """ Run the simulation for a given time
        """
        time *= self.parameters['time_conversion']

        cdef:
            long external_n = long(time/self.parameters['dt']/internal_steps)
            long internal_n = internal_steps
            float[:,::1] moments = np.zeros((external_n, 3), dtype=np.float32)

        times = self.times(moments, internal_steps)

        for i in range(external_n):
            for j in range(internal_n):
                self.previous = self.current
                self.step_func(self)

            moments[i][0] = self.current.m.x;
            moments[i][1] = self.current.m.y;
            moments[i][2] = self.current.m.z;

        return times, np.asarray(moments)


    def relax(self, steps=1000):
        """ Run the simulation until no torque is present above a threshold
        """
        cdef:
            float g0 = 0.0
            float g1 = 0.0
            long i
            long n = steps

        for i in range(n):
            self.previous = self.current
            self.step_func(self)

        g0 = self.energy(self.current.t, self.current.m)

        for i in range(n):
            self.previous = self.current
            self.step_func(self)

        g1 = self.energy(self.current.t, self.current.m)

        # Minimize the energy
        while g0 > g1:
            g0 = g1
            for i in range(n):
                self.previous = self.current
                self.step_func(self)
            g1 = self.energy(self.current.t, self.current.m)




    def times(self, moments, internal_steps=250):
        """ Returns an array of times that correspond to the moments array
        that is about to be run
        """
        return (self.t + internal_steps*self.parameters['dt']*(
                1 + np.arange(moments.shape[0])))/(
                self.parameters['time_conversion'])


    cdef float3 field(self, float t, float3 m):
        cdef float3 field
        return field


    cdef float3 torque(self, float t, float3 m):
        cdef float3 torque
        return torque


    cdef float energy(self, float t, float3 m):
        return energy.zeeman(m, self.field(t, m))


    def energy_surface(self, points=100, time=0.0):
        """ Returns a numpy array of vectors that have a length of the 
        energy at their orientation
        """
        cdef:
            long n = points
            float t = time
            float[::1] theta = np.linspace(0, np.pi, num=n, dtype=np.float32)
            float[::1] phi = np.linspace(-np.pi, np.pi, num=n, dtype=np.float32)
            float[:,::1] energies = np.zeros((points**2, 3), dtype=np.float32)
            float3 m
            long i, j, idx

        for i in range(n):
            for j in range(n):
                m.x = sin(theta[i])*cos(phi[j])
                m.y = sin(theta[i])*sin(phi[j])
                m.z = cos(theta[i])
                g = self.energy(t, m)
                idx = n*i + j
                energies[idx][0] = g*m.x
                energies[idx][1] = g*m.y
                energies[idx][2] = g*m.z

        return np.asarray(energies)




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




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Basic Kernel
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


cdef class BasicKernel(Kernel):


    def _load(self):
        """ Loads the parameters from the kernel (extending Kernel._load)
        """
        super(BasicKernel, self)._load()
        self.alpha = self.parameters['damping']
        self.hext = make_float3(self.parameters['hext'])
        self.Nd = make_float3(self.parameters['Nd'])


    cdef float3 field(self, float t, float3 m):
        return self.hext + field.demagnetization(m, self.Nd)


    cdef float3 torque(self, float t, float3 m):
        cdef float3 heff = self.field(t, m)
        return torque.landau_lifshitz(m, heff, self.alpha)


    property hext:
        def __get__(self): return make_array(self.hext)
        def __set__(self, float[::1] l): self.hext = make_float3(l)


    property Nd:
        def __get__(self): return make_array(self.Nd)
        def __set__(self, float[::1] l): self.Nd = make_float3(l)


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Anisotropy Kernel
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


cdef class AnisotropyKernel(BasicKernel):


    def _load(self):
        """ Loads the parameters from the kernel (extending BasicKernel._load)
        """
        super(AnisotropyKernel, self)._load()
        self.u = make_float3(self.parameters['u'])
        self.hu1 = self.parameters['hu1']
        self.hu2 = self.parameters['hu2']
        self.c1 = make_float3(self.parameters['c1'])
        self.c2 = make_float3(self.parameters['c2'])
        self.c3 = self.c1.cross(self.c2)
        self.hc1 = self.parameters['hc1']
        self.hc2 = self.parameters['hc2']


    cdef float3 field(self, float t, float3 m):
        cdef float3 heff = self.hext + field.demagnetization(m, self.Nd)
        heff += field.uniaxial_anisotropy(m, self.u, self.hu1, self.hu2)
        heff += field.cubic_anisotropy(m, self.c1, self.c2, self.c3, 
                    self.hc1, self.hc2)
        return heff


    cdef float3 torque(self, float t, float3 m):
        cdef float3 heff = self.field(t, m)
        return torque.landau_lifshitz(m, heff, self.alpha)


    property hext:
        def __get__(self): return make_array(self.hext)
        def __set__(self, float[::1] l): self.hext = make_float3(l)


    property Nd:
        def __get__(self): return make_array(self.Nd)
        def __set__(self, float[::1] l): self.Nd = make_float3(l)