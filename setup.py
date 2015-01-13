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
    author_email='clj72@cornell.edu',
    packages=['macrospin'],
    scripts=[],
    url='http://ralphgroup.lassp.cornell.edu/',
    license='LICENSE.txt',
    description='Macrospin simulations using Cython and CUDA for Python',
    long_description=open('README.md').read(),
    ext_modules=cythonize(extensions),
    package_data={'macrospin': ['*.pxd', '*.h']},
)