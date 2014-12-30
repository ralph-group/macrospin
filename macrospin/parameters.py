#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Parameter classes for dealing with units and normalization
#
# macrospin Python package
# Authors: Colin Jermain
# Copyright: 2014 Cornell University
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
class Parameters(object):
    pass

class CgsParameters(Parameters):

    def normalize(self):
        """ Returns a normalized set of the simulation parameters """
        pass


class MksParameters(Parameters):
    pass