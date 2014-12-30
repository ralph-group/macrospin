import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.gpuarray as gpua
from pycuda.compiler import SourceModule

try:
	import pycuda.autoinit
	import pycuda.driver as cuda
	import pycuda.gpuarray as gpua
	from pycuda.compiler import SourceModule
except ImportError:
	raise ImportError("This part of the macrospin package requires pycuda")
from simulation import Kernel
from pkg_resources import resource_string
import numpy as np
import constants as const 

def default_template():
	return resource_string(__name__, 'kernel.cu')

class CUDAKernel(Kernel):

	def __init__(self, Ms, damping, threads, blocks, timestep=5e-13):
		super(CUDAKernel, self).__init__(Ms, damping, timestep, 
			template=default_template())
		self.threads, self.blocks = threads, blocks
		self.moments = np.zeros((self.threads*self.blocks, 4), np.float32)

		self.context['definitions'] = """
		const int i = blockIdx.x*blockDim.x + threadIdx.x;

		// Moment being operated on
		float4 mloc = m[i];
		float4 heff = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
		"""

	def initial_moments(self, initial_moment):
		""" Sets all moments to a particular direction
		"""
		self.moments[:,:] = initial_moment + [0.0]

	def index_operations(self, index_operations):
		self._index_operations = index_operations

	def static_field(self, Hx, Hy, Hz):
		self.context['static_field'] = """
		heff = heff + make_float4({Hx}f, {Hy}f, {Hz}f, 0.0f);
		""".format(Hx=np.float32(Hx), Hy=np.float32(Hy), Hz=np.float32(Hz))

	def demagnetization(self, Nx, Ny, Nz):
		self.context['demagnetization'] = """
		const float4 N = make_float4({Nx}f, {Ny}f, {Nz}f, 0.0f);
  		heff = heff - mult(mloc, N);
		""".format(Nx=np.float32(Nx), Ny=np.float32(Ny), Nz=np.float32(Nz))

	def uniaxial_anisotropy(self, Hku):
		self.context['uniaxial_anisotropy'] = """
		heff.z = heff.z + {{ Hku }}*mloc.z;
		""".format(Hku=np.float32(Hku))

	def cubic_anisotropy(c1, c2, Hkc1, Hkc2=0.0):
		self.context['cubic_anisotropy'] = """
		float4 m2 = make_float4(mloc.x*mloc.x, mloc.y*mloc.y, mloc.z*mloc.z, 0.0f);
		  
		float hcx = {Hkc1}*(m2.y+m2.z) + {Hkc2}*(m2.y*m2.z);
		float hcy = {Hkc1}*(m2.x+m2.z) + {Hkc2}*(m2.x*m2.z);
		float hcz = {Hkc1}*(m2.x+m2.y) + {Hkc2}*(m2.x*m2.y);
		heff = heff - make_float4(-mloc.x*hcx, -mloc.y*hcy, -mloc.z*hcz, 0.0f);
		""".format(Hkc1=np.float32(Hkc1), Hkc2=np.float32(Hkc2))

	def spin_torque(self, Jc, Ms, thickness, polarization):
		Jcx, Jcy, Jcz = Jc
		torque_factor = ((const.hbar*0.1)/
			(2.0*const.ech*Ms**2*thickness*self.damping))*polarization
		self.context['spin_torque'] = """
		const float4 stt = {torque_factor}*Jc;
		float4 pxm = cross(stt,  mloc);
		hxm = hxm + pxm
		""".format(torque_factor=np.float32(-torque_factor))
		self.context['definitions'] += """
		float4 Jc = make_float4({Jcx}f, {Jcy}f, {Jcz}f, 0.0f);
		""".format(Jcx=np.float32(Jcx), Jcy=np.float32(Jcy), Jcz=np.float32(Jcz))

	def compile(self):
		self.gpu_moments = gpua.to_gpu(self.moments.astype(np.float32))

		self.model = SourceModule(self.render())
		self.evolve = self.model.get_function("evolve")
		self.normalize = self.model.get_function("normalize")

		self._i = 0

	def generator(self):
		while self._i > -1:
			for j in range(250):
				self.evolve(self.gpu_moments, 
							grid=(self.blocks,1), block=(self.threads,1,1))
			self.normalize(self.gpu_moments, grid=(self.blocks,1), block=(self.threads,1,1))
			if self._i % 10 == 0:
				self.moments = self.gpu_moments.get()
				self.callback(self.moments, self.timestep*self._i*250.0)
			self._i += 1
			yield