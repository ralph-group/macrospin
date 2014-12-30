#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Simulation class for threading Kernel objects
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
