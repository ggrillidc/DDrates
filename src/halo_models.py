#######################
#                     #
#     Halo models     #
#                     #
#######################



################
# Explanations #
################

#This program contains the Standard Halo Model


##########
# Import #
##########

#This part of the code imports the necessary Python libraries.

#Standard libraries
import math  as math
import numpy as np
from numpy import sqrt, sin, cos, pi, exp, heaviside, minimum, maximum
import scipy as scipy
import numericalunits as nu
from scipy.special import erf

nu.reset_units()

def SpeedDistIso(v,v0,v_esc):
    v_pec   = np.array([11.1,12.2,7.3]) * nu.km / nu.s #solar peculiar velocity
    v_LSR   = np.array([0.,220.,0.])  * nu.km / nu.s #local standard at rest velocity
    v_earth = np.sum((v_pec+v_LSR)**2)**0.5
    v_max   = v_esc + v_earth

    #Normalization constant
    v_ratio = v_esc / v0

    Nesc = erf(v_ratio) - 2 / sqrt(pi) * v_ratio * exp(-v_ratio**2)

    xmax = minimum(1, (v_esc**2 - v_earth**2 - v**2)/(2 * v_earth * v))

    etav = (Nesc * v / (pi**0.5 * v0 * v_earth) * (exp(-((v-v_earth)/v0)**2)
            - exp(-(v**2 + v_earth**2 + 2 * v * v_earth * xmax)/v0**2)))

    if v > v_max:
        etav = 0

    return etav

#Maximum observable dark matter velocity on Earth.
def v_max(v_esc):
    # solar peculiar velocity
    v_pec = np.array([11.1,12.2,7.3]) * nu.km / nu.s
    #local standard at rest velocity
    v_LSR = np.array([0.,220.,0.])  * nu.km / nu.s
    v_earth   = np.sum((v_pec+v_LSR)**2)**0.5
    v_max_tmp = v_esc + v_earth
    return v_max_tmp


class HaloModels:

    def __init__(self):
        self.v0 = 220  * nu.km / nu.s
        self.v_esc = 544  * nu.km / nu.s
        self.rho_DM = 0.3 * nu.GeV / nu.c0**2 / nu.cm**3

    def velocity_dist_SHM(self,v):
        return SpeedDistIso(v,self.v0,self.v_esc)


#end
