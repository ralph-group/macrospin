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


	def verify(self):
		""" Verifies that the required parameters are present
		"""
		raise NotImplementedError("Parameter subclass must define the verify method")


	def normalize(self):
		""" Normalizes the fields and times into those used for simulation,
		and returns modified Parameter object
		"""
		raise NotImplementedError("Parameter subclass must define the normalize method")


class CgsParameters(Parameters):


	def normalize_H(self, H):
		return np.asarray(H, dtype=np.float32)/self['Ms']

	def normalize(self):
		if 'gyromagnetic_ratio' not in self:
			self['gyromagnetic_ratio'] = constants.gyromagnetic_ratio

		# Rescale gyromagnetic ratio to Landau-Lifshitz version
		self['gyromagnetic_ratio'] = self['gyromagnetic_ratio']/(1.0 + self['damping']**2.)

		# Rescale time in terms of (gamma Ms)^-1
		self['time_conversion'] = self['gyromagnetic_ratio']*self['Ms']
		self['dt'] = self['dt']*self['time_conversion']

		# Rescale fields in terms of Ms
		for key, value in self.iteritems():
			if key.startswith("H"): # Assume it is a field
				self[key] = self.normalize_H(value)

		# Ensure initial moment is a numpy array
		self['m0'] = np.asarray(self['m0'], dtype=np.float32)

		self['Nd'] = np.asarray(self['Nd'], dtype=np.float32)

		return self


class MksParameters(Parameters):
	pass