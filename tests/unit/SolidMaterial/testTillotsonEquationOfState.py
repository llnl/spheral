#ATS:test(SELF, label="Tillotson EOS unit tests.")
# Unit tests of the Tillotson equation of state.  This is just checking the coding by verifying results
# against direct equations spit out by Mathematica.
import unittest
from math import *
from SolidSpheral1d import *

#===============================================================================
# Unit tests.
#===============================================================================
class TestTillotsonEquationOfState(unittest.TestCase):

    #===========================================================================
    # setUp
    #===========================================================================
    def setUp(self):
        self.nsample = 100
        self.Ptol = 1.0e-10

        self.logRhoMin = log(1e-5)
        self.logRhoMax = log(1e5)
        self.drho = (self.logRhoMax - self.logRhoMin)/self.nsample

        self.logEpsMin = log(1e-10)
        self.logEpsMax = log(1e10)
        self.deps = (self.logEpsMax - self.logEpsMin)/self.nsample

        # I largely lift these parameters from Saito et al. (2008) for
        # the Asteroid dolomite material, converted here for MKS units.
        self.rho0 = 3.0e3
        self.etamin, self.etamax = 1e-10, 1e10
        self.a = 0.5
        self.b = 0.6
        self.A = 91.1e9
        self.B = 32.1e9
        self.alpha = 5.0
        self.beta = 5.0
        self.eps0 = 10.0e6
        self.epsLiquid = 250.0e6
        self.epsVapor = 1.4e9
        self.atomicWeight = 20.12  # Taken from Murty 1962 for Limestone.
     
        self.eos = TillotsonEquationOfStateMKS(self.rho0,
                                               self.etamin,
                                               self.etamax,
                                               self.etamin,
                                               self.etamax,
                                               self.a,
                                               self.b,
                                               self.A,
                                               self.B,
                                               self.alpha,
                                               self.beta,
                                               self.eps0,
                                               self.epsLiquid,
                                               self.epsVapor,
                                               self.atomicWeight)

        self.nodes = makeFluidNodeList("test nodes", self.eos, numInternal=1)

    #===========================================================================
    # tearDown
    #===========================================================================
    def tearDown(self):
        del self.nodes, self.eos

    #===========================================================================
    # Pressure analytic answer.
    #===========================================================================
    def P1(self, rhoi, epsi, etai, mui):
        return (self.a + self.b/(1.0 + epsi/(self.eps0*etai*etai)))*rhoi*epsi + self.A*mui + self.B*mui*mui

    def P2(self, rhoi, epsi, etai, mui):
        return (self.a + self.b/(1.0 + epsi/(self.eps0*etai*etai)))*rhoi*epsi + self.A*mui

    def P3(self, rhoi, epsi, etai, mui):
        p2 = self.P2(rhoi, self.epsLiquid, etai, mui)
        p4 = self.P4(rhoi, self.epsVapor, etai, mui)
        return p2 + (p4 - p2)*(epsi - self.epsLiquid)/(self.epsVapor - self.epsLiquid)

    def P4(self, rhoi, epsi, etai, mui):
        return (self.a*rhoi*epsi +
                (self.b*rhoi*epsi/(1.0 + epsi/(self.eps0*etai*etai)) +
                 self.A*mui*exp(self.beta*(1.0 - 1.0/etai)))*exp(-self.alpha*(1.0 - 1.0/etai)**2))

    def Pans(self, rhoi, epsi):
        etai = self.eos.boundedEta(rhoi)
        mui = etai - 1.0
        rho = etai*self.rho0
        phi = self.b/(1.0 + epsi/(self.eps0*etai*etai));
        chi = 1.0/etai - 1.0;
        if mui >= 0.0:
            return (self.a + phi)*rho*epsi + self.A*mui + self.B*mui*mui
        elif epsi <= self.epsLiquid:
            return (self.a + phi)*rho*epsi + self.A*mui # + self.B*mui*mui
        elif epsi >= self.epsVapor:
            return self.a*rho*epsi + (phi*rho*epsi + self.A*mui*exp(-self.beta*chi))*exp(-self.alpha*chi*chi)
        else:
            P2 = (self.a + phi)*rho*epsi + self.A*mui  #  + self.B*mui*mui
            P4 = self.a*rho*epsi + (phi*rho*epsi + self.A*mui*exp(-self.beta*chi))*exp(-self.alpha*chi*chi)
            return P2 + (P4 - P2)*(epsi - self.epsLiquid)/(self.epsVapor - self.epsLiquid)

    #===========================================================================
    # dPdrho analytic answer.
    #===========================================================================
    def dPdrho1(self, rhoi, epsi, etai, mui):
        return (self.dPdrho2(rhoi, epsi, etai, mui) +
                2.0*self.B*(rhoi - self.rho0)/(self.rho0**2))

    def dPdrho2(self, rhoi, epsi, etai, mui):
        return (self.a*epsi +
                self.A/self.rho0 +
                2.0*self.b*epsi**2*self.eps0*(rhoi*self.rho0)**2/(self.eps0*rhoi**2 + epsi*self.rho0**2)**2 +
                self.b*epsi*self.eps0*rhoi**2/(self.eps0*rhoi**2 + epsi*self.rho0**2))

    def dPdrho3(self, rhoi, epsi, etai, mui):
        dp2dhro = self.dP2drho(rhoi, epsi, etai, mui)
        dp4drho = self.dP4drho(rhoi, epsi, etai, mui)
        return dp2drho + (dp4drho - dp2drho)*(epsi - self.epsLiquid)/(self.epsVapor - self.epsLiquid)

    def dPdrho4(self, rhoi, epsi, etai, mui):
        return (self.a*epsi -
                (2*self.alpha*exp(-self.alpha*((rhoi - self.rho0)/rhoi)**2) * self.rho0*(1.0 - 1.0/etai) *
                 (self.A*exp(self.beta*(1.0 - 1.0/etai))*(1.0/etai - 1.0) +
                  self.b*epsi*self.eps0*rhoi**3/(self.eps0*rhoi**2 + epsi*self.rho0**2)))/self.rho**2 +
                exp(-self.alpha*((rhoi - self.rho0)/rhoi)**2)*
                (self.A*exp(self.beta*(1.0 - 1.0/etai))*(rhoi**2 + self.beta*rhoi*self.rho0 - self.beta*self.rho0**2)/
                 (rhoi**2*self.rho0) +
                 self.b*epsi*self.eps0*rhoi**2*(self.eps0*rhoi**2 + 3.0*eps*self.rho0**2)/
                 (self.eps0*rhoi**2 + epsi*self.rho0**2)**2))

    def dPdrhoans(self, rhoi, epsi):
        etai = self.eos.boundedEta(rhoi)
        mui = etai - 1.0
        if mui >= 0.0:
            print("Region 1")
            return self.P1(rhoi, epsi, etai, mui)
        elif epsi <= self.epsLiquid:
            print("Region 2")
            return self.P2(rhoi, epsi, etai, mui)
        elif epsi <= self.epsVapor:
            print("Region 3")
            return self.P3(rhoi, epsi, etai, mui)
        else:
            print("Region 4")
            return self.P4(rhoi, epsi, etai, mui)

    #===========================================================================
    # a
    #===========================================================================
    def testa(self):
        assert self.eos.a == self.a

    #===========================================================================
    # b
    #===========================================================================
    def testb(self):
        assert self.eos.b == self.b

    #===========================================================================
    # A
    #===========================================================================
    def testA(self):
        assert self.eos.A == self.A

    #===========================================================================
    # B
    #===========================================================================
    def testB(self):
        assert self.eos.B == self.B

    #===========================================================================
    # alpha
    #===========================================================================
    def testalpha(self):
        assert self.eos.alpha == self.alpha

    #===========================================================================
    # beta
    #===========================================================================
    def testbeta(self):
        assert self.eos.beta == self.beta

    #===========================================================================
    # eps0
    #===========================================================================
    def testeps0(self):
        assert self.eos.eps0 == self.eps0

    #===========================================================================
    # epsLiquid
    #===========================================================================
    def testepsLiquid(self):
        assert self.eos.epsLiquid == self.epsLiquid

    #===========================================================================
    # epsVapor
    #===========================================================================
    def testepsVapor(self):
        assert self.eos.epsVapor == self.epsVapor

    #===========================================================================
    # atomicWeight
    #===========================================================================
    def testatomicWeight(self):
        assert self.eos.atomicWeight == self.atomicWeight

    #===========================================================================
    # rho
    #===========================================================================
    def rho(self, i):
        return exp(self.logRhoMin + (i + 0.5)*self.drho)

    #===========================================================================
    # eps
    #===========================================================================
    def eps(self, i):
        return exp(self.logEpsMin + (i + 0.5)*self.deps)

    #===========================================================================
    # pressure
    #===========================================================================
    def testPressure(self):
        rhof = ScalarField("rho", self.nodes)
        epsf = ScalarField("eps", self.nodes)
        Pf = ScalarField("pressure", self.nodes)
        for irho in range(self.nsample):
            rhof[0] = self.rho(irho)
            for ieps in range(self.nsample):
                epsf[0] = self.eps(ieps)
                self.eos.setPressure(Pf, rhof, epsf)
                Pi = Pf[0]
                P0 = self.Pans(rhof[0], epsf[0])
                eta = self.eos.boundedEta(rhof[0])
                mu = eta - 1.0
                phi = self.eos.computePhi(eta, epsf[0])
                P2 = self.eos.computeP2(phi, mu, rhof[0], epsf[0])
                self.assertTrue(fuzzyEqual(Pi, P0, self.Ptol),
                                "Pressure do not match:  P(%g, %g) = %g != %g\n P1=(%g,%g) P2=(%g,%g), P4=(%g,%g)\n eta=%g mu=%g phi=%g" % 
                                (rhof[0], epsf[0], Pi, P0,
                                 self.eos.computeP1(mu, P2), self.P1(rhof[0], epsf[0], eta, mu),
                                 P2, self.P2(rhof[0], epsf[0], eta, mu),
                                 self.eos.computeP4(phi, mu, eta, rhof[0], epsf[0]), self.P4(rhof[0], epsf[0], eta, mu),
                                 eta, mu, phi))
        return

    #===========================================================================
    # dPdrho
    #===========================================================================
    # def testdPdrho(self):
    #     for irho in xrange(self.nsample):
    #         rhoi = self.rho(irho)
    #         for ieps in xrange(self.nsample):
    #             epsi = self.eps(ieps)
    #             dPdrhoi = self.eos.computeDPDrho(rhoi, epsi)
    #             dPdrho0 = self.dPdrhoans(rhoi, epsi)
    #             self.assertTrue(fuzzyEqual(dPdrhoi, dPdrho0, self.Ptol),
    #                             "dP/drho does not match:  dP/drho(%g, %g) = %g != %g" % (rhoi, epsi, dPdrhoi, dPdrho0))
    #     return

#===============================================================================
# Run the suckers.
#===============================================================================
if __name__ == "__main__":
    unittest.main()
