#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Simulation class for threading Kernel objects
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
import numpy as np
from threading import Thread

class Simulation(object):
    """ Simulation object coordinates the Thread that runs the kernel
    """
    FAILED, STOPPED, QUEUED, RUNNING = 0, 1, 2, 3

    def __init__(self, kernel):
        self.kernel = kernel
        self.status = Simulation.QUEUED
        self.thread = None

    def run(self, time, timeout=None, thread_timeout=1e5):
        """ Runs the simulation for a given simulation time in seconds
        """
        if self.status == Simulation.RUNNING:
            raise Exception("Simulation is already running")
        elif self.status == Simulation.FAILED:
            raise Exception("Can not run a simulation that has failed")
        self.status = Simulation.RUNNING
        self.thread = Thread(target=self.kernel, kwargs={
                            'time': time, 'timeout': timeout})
        self.thread.start()

    def isRunning(self):
        """ Returns True if an existing Thread is running
        """
        if self.thread is None:
            return False
        else:
            return self.thread.isAlive()

    def wait(self, thread_timeout):
        self.thread.join(thread_timeout)
        self.status = Simulation.STOPPED

    def stop(self):
        """ Stops the simulation even if running
        """
        self.kernel.stop()
        
    def resume(self, time, thread_timeout=1e5):
        """ Resumes the simulation from the last moment orientation
        """
        self.run(time, None, thread_timeout)

    def stabilize(self, timeout=1e-4, thread_timeout=1e5):
        """ Runs the simulation until the moment is stable or the timeout 
        is reached in simulation time
        """
        self.run(None, timeout, thread_timeout)


class FieldSweep(object):


    def __init__(self, kernel, fields):
        self.kernel = kernel
        self.fields = fields


    @staticmethod
    def loop(kernel, direction='x', start_field=-1e3, end_field=1e3, 
             points=1e3, reverse=True):
        """ Returns a FieldSweep object with fields that go from the start_field
        to the end_field along a specific direction, with the default option to
        also include the reverse
        """
        H = np.linspace(start_field, end_field, num=points, dtype=np.float32)
        coordinates = {'x': 0, 'y': 1, 'z': 2}

        if direction not in coordinates:
            raise ValueError("Field sweep direction must be either x, y, or z")

        if reverse:
            fields = np.zeros((2*points, 3), dtype=np.float32)
            fields[:points,coordinates[direction]] = H
            fields[points:,coordinates[direction]] = H[::-1] # Reversed view
        else:
            fields = np.zeros((points, 3), dtype=np.float32)
            fields[:,coordinates[direction]] = H

        return FieldSweep(kernel, fields)


    def run(self, threshold=1e-3):
        """ Runs through each field and stabilizes the moment, returning
        the fields, stabilization time, and moment orientation
        """
        size = self.fields.shape[0]
        times = np.zeros((size, 1), dtype=np.float32)
        moments = np.zeros((size, 3), dtype=np.float32)

        self.kernel.reset()

        for i, field in enumerate(self.fields):
            ti = self.kernel.t_sec
            self.kernel.parameters['Hext'] = self.kernel.parameters.normalize_H(field)
            self.kernel.stabilize(threshold)
            times[i] = self.kernel.t_sec - ti
            moments[i] = self.kernel.m

        return self.fields, times, moments

