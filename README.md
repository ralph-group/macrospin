# macrospin Python package #

The `macrospin` package simulates isolated macrospins in a simple and efficient way. The Landau-Lifshitz equation is solved by either a C++-like Cython kernel or a CUDA kernel. Both are accessible directly from Python to allow easy simulations in and out of the IPython Notebook.

## Basic example of macrospin in applied field ##

```python
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
```