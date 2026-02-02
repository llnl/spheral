#-------------------------------------------------------------------------------
# ShadowPalphaPorosity
#
# Provides convenient constructors for the P-alpha porosity model
#-------------------------------------------------------------------------------
import types
import scipy.integrate
from spheralDimensions import spheralDimensions
from CaptureStdout import helpString
from MaterialPropertiesLib import SpheralMaterialPropertiesLib
from SpheralCompiledPackages import *

dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Factory function to generate our Python enhanced PalphaPorosity model,
# including using scipy's fsolve to find alphae as needed.
#-------------------------------------------------------------------------------
def _PalphaPorosityFactory(ndim):
    CXXPalphaPorosity = eval(f"PalphaPorosity{ndim}d")
    ScalarField = eval(f"ScalarField{ndim}d")

    # The class we will return in place of the plain C++ one
    class PalphaPorosity(CXXPalphaPorosity):

        #-----------------------------------------------------------------------
        # Constructor
        #-----------------------------------------------------------------------
        def __init__(self,
                     nodeList,
                     phi0,
                     Pe,
                     Pt,
                     Ps,
                     alphat,
                     n1,
                     n2,
                     cS0,
                     c0,
                     rhoS0 = None,
                     alphae = None,
                     jutziStateUpdate = True):
            """Construct a P-alpha porosity model.  Valid arguments are

  *         nodeList: the solid NodeList whose nodes this porosity model will apply to                                           
  *             phi0: (scalar or Field) initial porosity                                                                         
  *               Pe: (scalar) Elastic pressure threshold                                                                        
  *               Pt: (scalar) Transition pressure (Pe <= Pt)                                                                    
  *               Ps: (scalar) Solid transition pressure (from fit, Pt <= Ps)                                                    
  *           alphat: (scalar) Transition distension                                                                             
  *               n1: (scalar) Fitted exponent for plastic distention evolution                                                  
  *               n2: (scalar) Fitted exponent for plastic distention evolution                                                  
  *              cS0: (scalar) Reference sound speed at full density                                                             
  *               c0: (scalar or Field) Reference sound speed at initial porosity                                                
  *            rhoS0: (scalar or None) Reference (solid) density.  If None, then reference value looked up from equation of state
  *           alphae: (scalar or None) Elastic transition distension. If None, then iteratively solved for.
  * jutziStateUpdate: (bool, default True) Update state (deviatoric stress and damage) as described in Jutzi 2008
"""

            if rhoS0 is None:
                try:
                    rhoS0 = nodeList.equationOfState().referenceDensity
                except:
                    raise RuntimeError("Unable to extract reference density for PalphaPorosity")

            # Find the starting maximum porosity (maximum distension) and corresponding minimum initial sound speed
            if isinstance(phi0, ScalarField):
                alpha0 = 1.0/(1.0 - phi0.max())
                c0min = c0.min()
            else:
                alpha0 = 1.0/(1.0 - phi0)
                c0min = c0

            # Now find alphae = alpha(Pe)
            if alphae is None:
                if Pe > 0.0:
                    K0 = rhoS0*cS0*cS0
                    alphae = alpha0   # Starting point
                    last_alphae = 0.0
                    iter = 0
                    def DalphaDP_elastic(P, alpha):
                        h = 1.0 + (alpha - 1.0)*(c0min - cS0)/max(1.0e-20, cS0*(alphae - 1.0))
                        return alpha*alpha/K0*(1.0 - 1.0/(h*h))
                    while abs(alphae - last_alphae) > 1.0e-10 and iter < 1000:
                        iter += 1
                        last_alphae = alphae
                        #print(" --> ", iter, alphae, c0min, cS0, K0)
                        alphae = scipy.integrate.solve_ivp(DalphaDP_elastic,
                                                           #args = (c0min, cS0, K0, alphae),
                                                           t_span = [0.0, Pe],
                                                           y0 = [alpha0],
                                                           t_eval = [Pe]).y[0][0]
                    #print("alphae: ", alphae, last_alphae, iter)
                else:
                    alphae = alpha0

            if alphat is None:
                alphat = alphae    # Reduces the plastic crush curve to Eq. 8 from Jutzi 2008

            # Remember the alpha0 max and c0min for use in the convenience crushCurve method
            self._alpha0max = alpha0
            self._c0min = c0min

            # Now build the actual C++ object
            CXXPalphaPorosity.__init__(self,
                                       nodeList = nodeList,
                                       phi0 = phi0,
                                       Pe = Pe,
                                       Pt = Pt,
                                       Ps = Ps,
                                       alphae = alphae,
                                       alphat = alphat,
                                       n1 = n1,
                                       n2 = n2,
                                       cS0 = cS0,
                                       c0 = c0,
                                       rhoS0 = rhoS0,
                                       jutziStateUpdate = jutziStateUpdate)
            return
    
        #-----------------------------------------------------------------------
        # Compute the crush-curve on demand as a function of pressure
        #-----------------------------------------------------------------------
        def crushCurve(self, P):
            if (P >= self.Ps):

                # Solid limit, all porosity squeezed out
                return 1.0

            elif P < self.Pe:
                assert self._c0min < self.cS0

                # Elastic regime -- integrate elastic Dalpha/DP equation
                stuff = scipy.integrate.solve_ivp(self.Dalpha_elasticDP,
                                                  args = (self._c0min, self.cS0, self.K0, self.alphae),
                                                  t_span = [0.0, P],
                                                  y0 = [self._alpha0max],
                                                  t_eval = [P])
                return stuff.y[0][0]

            else:

                # Plastic limit
                if P < self.Pt:
                    return (self.alphae - self.alphat) * ((self.Pt - P)/(self.Pt - self.Pe))**self.n1 + (self.alphat - 1.0) * ((self.Ps - P)/(self.Ps - self.Pe))**self.n2 + 1.0
                else:
                    return (self.alphat - 1.0) * ((self.Ps - P)/(self.Ps - self.Pe))**self.n2 + 1.0

        #-----------------------------------------------------------------------
        # Reproduce the Dalpha/DP relationship in the elastic regime
        #-----------------------------------------------------------------------
        def DalphaDP_elastic(self, P, alpha,
                             c0min, cS0, K0, alphae):
            h = 1.0 + (alpha - 1.0)*(c0min - cS0)/(cS0*(alphae - 1.0))
            return alpha*alpha/K0*(1.0 - 1.0/(h*h))

    return PalphaPorosity

#-------------------------------------------------------------------------------
# Create the dimension specific factories.  These are the ones you actually use.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec(f"PalphaPorosity{ndim}d = _PalphaPorosityFactory({ndim})")
