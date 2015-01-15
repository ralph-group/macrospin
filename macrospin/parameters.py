#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter classes for dealing with units and normalization
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014-2015 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
from macrospin import constants
import numpy as np

class Parameters(dict):
	pass


class CgsParameters(Parameters):
	

	def normalize_field(self, field_name):
		if isinstance(field_name, str):
			value = self[field_name]
		else:
			value = field_name
		return np.asarray(value, dtype=np.float32)/self['Ms']


	def normalize_energy(self, energy_name):
		if isinstance(energy_name, str):
			value = self[energy_name]
		else:
			value = energy_name
		return 2.0*self[energy_name]/self['Ms']**2.0


class MksParameters(Parameters):
	pass


class NormalizedParameters(Parameters):


	def __init__(self, parameters):

		self['damping'] = parameters['damping']
		self['m0'] = np.asarray(parameters['m0'], dtype=np.float32)
		self['Ms'] = parameters['Ms']
		self['dt'] = parameters['dt']

		if 'u' in parameters:
			self['u'] = np.asarray(parameters['u'], dtype=np.float32)
		else:
			self['c2'] = np.zeros((3,), dtype=np.float32)
		if 'c1' in parameters:
			self['c1'] = np.asarray(parameters['c1'], dtype=np.float32)
		else:
			self['c1'] = np.zeros((3,), dtype=np.float32)
		if 'c2' in parameters:
			self['c2'] = np.asarray(parameters['c2'], dtype=np.float32)
		else:
			self['c2'] = np.zeros((3,), dtype=np.float32)

		if 'gyromagnetic_ratio' not in parameters:
			self['gyromagnetic_ratio'] = constants.gyromagnetic_ratio

		# Rescale gyromagnetic ratio to Landau-Lifshitz version
		self['gyromagnetic_ratio'] = self['gyromagnetic_ratio']/(1.0 + self['damping']**2.)

		# Rescale time in terms of (gamma Ms)^-1
		self['time_conversion'] = self['gyromagnetic_ratio']*self['Ms']
		self['dt'] = self['dt']*self['time_conversion']

		if 'Nd' in parameters:
			self['Nd'] = np.asarray(parameters['Nd'], dtype=np.float32)
		else:
			self['Nd'] = np.zeros((3,), dtype=np.float32)

		if 'Hext' in parameters: self['hext'] = parameters.normalize_field('Hext')

		if 'Ku1' in parameters: self['hu1'] = parameters.normalize_energy('Ku1')
		else: self['hu1'] = 0.0

		if 'Ku2' in parameters:	self['hu2'] = parameters.normalize_energy('Ku2')
		else: self['hu2'] = 0.0

		if 'Kc1' in parameters:	self['hc1'] = parameters.normalize_energy('Kc1')
		else: self['hc1'] = 0.0

		if 'Kc2' in parameters: self['hc2'] = parameters.normalize_energy('Kc2')
		else: self['hc2'] = 0.0
