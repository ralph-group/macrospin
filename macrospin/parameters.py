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
	pass


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

		if isinstance(parameters, CgsParameters): 
			self.normalize_cgs(parameters)
		elif isinstance(parameters, MksParameters): 
			self.normalize_mks(parameters)


	def normalize_cgs(self, parameters):
		""" Fills this dictionary with normalized parameters from CGS units
		"""
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

		if 'Hext' in parameters:
			self['hext'] = np.asarray(parameters['Hext'], dtype=np.float32)/self['Ms']

		if 'Ku1' in parameters:
			self['hu1'] = 2.0*parameters['Ku1']/self['Ms']**2.0
		else:
			self['hu1'] = 0.0
		if 'Ku2' in parameters:
			self['hu2'] = 2.0*parameters['Ku2']/self['Ms']**2.0
		else:
			self['hu2'] = 0.0
		if 'Kc1' in parameters:
			self['hc1'] = 2.0*parameters['Kc1']/self['Ms']**2.0
		else:
			self['hc1'] = 0.0
		if 'Kc2' in parameters:
			self['hc2'] = 2.0*parameters['Kc2']/self['Ms']**2.0
		else:
			self['hc2'] = 0.0


	def normalize_mks(self, parameters):
		""" Fills this dictionary with normalized parameters from MKS units
		"""
		pass