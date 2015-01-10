from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension("macrospin.kernels",
        ["macrospin/kernels.pyx", "macrospin/cpp/kernels.cpp", "macrospin/cpp/solvers.cpp"], 
        language='c++',
        include_dirs=[np.get_include(), 'macrospin/cpp/'],
        library_dirs=['macrospin/cpp/'],
        extra_compile_args=['-std=c++11'], # this line might break Windows/MacOSX 
    ),
]

setup(
    name='macrospin',
    version='0.1',
    author='Colin Jermain',
    author_email='clj72@cornell.edu',
    packages=['macrospin'],
    scripts=[],
    url='http://ralphgroup.lassp.cornell.edu/',
    license='LICENSE.txt',
    description='Macrospin simulations using Cython and CUDA for Python',
    long_description=open('README.md').read(),
    ext_modules=cythonize(extensions),

)