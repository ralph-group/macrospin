from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("macrospin.kernels",
        ["macrospin/kernels.pyx"], 
        language='c++',
        include_dirs=[np.get_include(), '.'],
    ),
]

setup(
    name='macrospin',
    version='0.1',
    author='Colin Jermain',
    packages=['macrospin'],
    scripts=[],
    url='http://ralphgroup.lassp.cornell.edu/',
    license='MIT License',
    description='Macrospin simulations using Cython and CUDA for Python',
    long_description=open('README.md').read(),
    ext_modules=cythonize(extensions),
    package_data={'macrospin': ['*.pxd', '*.h']},
)