from macrospin.parameters import CgsParameters
from macrospin.kernels import BasicKernel
from macrospin.simulation import Simulation

parameters = CgsParameters({
	'Ms': 140, # emu/cc
	'dt': 5e-13, # sec
	'damping', 0.01,
	'H_ext': [0., 1e3, 0.], # Oe
})

data_filename = "example.csv"

kernel = BasicKernel(parameters)

simulation = Simulation(kernel, data_filename)
simulation.run(1e-7)

simulation.wait()

