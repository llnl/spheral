#-------------------------------------------------------------------------------
# A mock physics package to mess around with the CRKSPH corrections.
#-------------------------------------------------------------------------------
from Spheral1d import *

class CRKSPH_mod_package(Physics):

    def __init__(self):
        Physics.__init__(self)
        return

    def evaluateDerivatives(self, t, dt, db, state, derivs):
        return

    def dt(self, db, state, derivs, t):
        return (1e100, "No vote")

    def registerState(self, dt, state):
        return

    def registerDerivatives(self, db, derivs):
        return

    def label(self):
        return "CRKSPH_mod_package"

    def initialize(self, t, dt, db, state, derivs):

        # Grab the CRKSPH arrays.
        A0_fl = state.scalarFields(HydroFieldNames.A0_CRKSPH)
        A_fl = state.scalarFields(HydroFieldNames.A_CRKSPH)
        B_fl = state.vectorFields(HydroFieldNames.B_CRKSPH)

        A0 = A0_fl[0]
        A = A_fl[0]
        B = B_fl[0]

        print("A", A.internalValues())
        return

