#########################
#                       #
#    Brehmstrahlung     #
#                       #
#########################



################
# Explanations #
################

#This program contains the DM rate exploiting the elastic nuclear recoil detected by bremsstrahlung,
#see [arxiv:1607.01789v2]

##########
# Import #
##########

#This part of the code imports the necessary Python libraries.

#Standard libraries
import math  as math
import numericalunits as nu
import numpy as np
from numpy import sqrt, sin, cos, pi, exp, heaviside, minimum, maximum
import scipy as scipy
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from scipy.integrate import quad

import halo_models as hm
import DM_NR as DMNR
import DM_LAr as LAr


# Load the X-ray form factor
def f_inter(fn):
    dataf  = np.loadtxt('../data/bremsstrahlung/' + fn + '.txt')
    x = dataf[:,0]
    y = dataf[:,1]
    return interp1d(x, y, fill_value='extrapolate')


f1_Ar = f_inter('f1_Ar')
f2_Ar = f_inter('f2_Ar')

#Minimum velocity in order to emit a photon with energy Eb from arxiv:1607.01789
def vmin_bs(Eb, mDM):
    return sqrt(2 * Eb / DMNR.mu_T(mDM))

#Brehmstrahlung boundaries of the recoil energy
#Eb: photon energy
#v: DM velocity
#mDM: DM mass

def ER_max(Eb, v, mDM):
    ret = (DMNR.mu_T(mDM)**2 * v**2 / DMNR.mT()
        * (1 - vmin_bs(Eb, mDM)**2 / (2 * v**2)
        + (1 - vmin_bs(Eb, mDM)**2 / v**2)**0.5))
    return ret

def ER_min(Eb, v, mDM):
    ret = (DMNR.mu_T(mDM)**2 * v**2 / DMNR.mT()
        * (1 - vmin_bs(Eb, mDM)**2 / (2 * v**2)
        - (1 - vmin_bs(Eb, mDM)**2 / v**2)**0.5))
    return ret


#Differential brehmstrahlung DM-nucleus cross section ($d^2\sigma/dE_R/d\omega$),
#from arxiv:1607.01789v2
def dsigmadERdEb_bs(Eb, ER, v, mDM, sigma_nucleon):
    halo_model = hm.HaloModels()
    v_max = hm.v_max(halo_model.v_esc)
    if Eb > DMNR.mu_T(mDM) * v_max**2 / 2.:
        ret = 0.
    else:
        # X-ray form factor
        F_atomic = np.abs(f1_Ar(Eb / nu.keV) + 1j * f2_Ar(Eb / nu.keV))
        ret = (4 * nu.alphaFS / (3 * np.pi * Eb) * ER / (DMNR.mT() * nu.c0**2)
            * F_atomic**2 * DMNR.dsigmadER(ER, v, mDM, sigma_nucleon))

    return ret

#Differential brehmstrahlung DM-nucleus cross section ($d^2\sigma/d\omega$),
#from arxiv:1607.01789v2
def dsigmadEb(Eb, v, mDM, sigma_nucleon):

    def dsdERdEb(ER):
        return dsigmadERdEb_bs(Eb, ER, v, mDM, sigma_nucleon)

    ret = quad(dsdERdEb, ER_min(Eb, v, mDM), ER_max(Eb, v, mDM))[0]

    return ret

#Differential rate per unit detector mass and recoil energy of the
#Bremsstrahlung process.
def rate_bremsstrahlung(Eb, mDM, sigma_nucleon):

    halo_model = hm.HaloModels()
    rhoDM = halo_model.rho_DM
    v_min = vmin_bs(Eb, mDM)
    v_max = hm.v_max(halo_model.v_esc)

    if v_min >= v_max:
        return 0

    def sigma(v):
        return (dsigmadEb(Eb, v, mDM, sigma_nucleon) * v * halo_model.velocity_dist_SHM(v))

    ret = rhoDM / mDM / DMNR.mT() * quad(sigma,v_min,v_max)[0]

    return ret


#Differential rate per unit detector mass and recoil energy of the
#Bremsstrahlung process.
def rate_bremsstrahlungN(Ne, mDM, sigma_nucleon):
    Eb         = LAr.Eer(Ne)
    halo_model = hm.HaloModels()
    rhoDM = halo_model.rho_DM
    v_min = vmin_bs(Eb, mDM)
    v_max = hm.v_max(halo_model.v_esc)

    if v_min >= v_max:
        return 0

    def sigma(v):
        return (dsigmadEb(Eb, v, mDM, sigma_nucleon) * v * halo_model.velocity_dist_SHM(v))

    rate = rhoDM / mDM / DMNR.mT() * quad(sigma,v_min,v_max)[0]
    rateN = LAr.dEerdN(Ne)*rate

    return rateN


def Nevents(mDM, sigma_nucleon):

    ret = quad(lambda x: rate_bremsstrahlungN(x, mDM, sigma_nucleon),7,50)[0]

    return ret * (nu.kg * nu.day)

#end
