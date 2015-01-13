# macrospin Python package #

The `macrospin` package simulates isolated macrospins in a simple and efficient way. The Landau-Lifshitz equation is solved by either a C++ kernel or a CUDA kernel. Both are accessible directly from Python to allow easy simulations in and out of the IPython Notebook.

## Basic example of macrospin in applied field ##

```python
from macrospin.parameters import CgsParameters
from macrospin.kernels import BasicKernel

parameters = CgsParameters({
    'Ms': 140, # Saturation Magnetization (emu/cc)
    'dt': 5e-13, # Timestep (sec)
    'damping': 0.01, # Gilbert damping
    'Hext': [0., 1e3, 0.], # External field (Oe)
    'm0': [-0.999, 0.001, 0.001], # Initial moment (normalized)
    'Nd': [0, 0, 0], # Demagnetization diagonal tensor elements
})

kernel = BasicKernel(parameters)

# times: Numpy array of simulation times
# moments: Numpy array of moment orientations
times, moments = kernel.run(time=1e-7)
```