##############################################
#                                            #
#     Nuclear recoil functions (only SI)     #
#                                            #
##############################################



################
# Explanations #
################

#This program contains the DM rate for Spin Independent scattering for Ar detectors

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
from scipy.integrate import quad

import halo_models as hm
import DM_LAr as LAr

#Mass of nucleus for Ar (target)
def mT():
    atomic_mass = 39.948
    return atomic_mass * nu.amu

#Reduced mass for a system with mass m1 and m2
def mu(m1, m2):
    return m1 * m2 / (m1 + m2)

#DM-nucleus reduced mass
def mu_T(mDM):
    return mu(mDM, mT())

#Minimum DM velocity for an elastic recoil ER and a DM mass mDM
def vmin_el(ER, mDM):
    y = np.sqrt(mT() * ER / (2 * mu_T(mDM)** 2))
    return y

#Spherical Bessel function of the first kind
def SphericalBesselJ(z):
    return sin(z) / z**2 - cos(z) / z

#Helm Form factor for Ar, correct with nu
def FHelm(ER):
    atomic_mass = 39.948
    fm = 5.06329 / nu.GeV
    mn = 0.938 * nu.GeV
    qER = sqrt(2 * mn * atomic_mass * ER)
    s = 0.9 * fm
    rn = 1.14 * atomic_mass**(1/3) * fm
    result = 3 * SphericalBesselJ(qER*rn) / (qER * rn) * exp(-(s * qER)**2 / 2)
    return result

#Maximum kinetic nuclear recoil energy
#mDM:        DM mass
#m_nucleus : nucleus mass
#v:          DM speed
def Emax(mDM, v, m_nucleus):
    return 2 * mu(mDM, m_nucleus)**2 * v**2 / m_nucleus

#Differential elastic DM-nucleus cross section (dependent on recoil energy and
#DM-earth speed v)
#ER:            recoil energy
#v:             DM speed (earth/detector frame)
#mDM:           mass of DM
#sigma_nucleon: DM-nucleon cross-section

def dsigmadER(ER, v, mDM, sigma_nucleon):
    atomic_mass = 39.948
    sigma_nucleus = (sigma_nucleon * (mu_T(mDM) / mu(nu.amu, mDM))**2 * atomic_mass**2)
    result = (sigma_nucleus / Emax(mDM, v, mT()) * FHelm(ER)**2)
    return result

#Differential rate per unit detector mass and recoil energy of elastic SI scattering
#ER:            recoil energy
#mDM:           mass of DM
#sigma_nucleon: DM-nucleon cross-section

def rate_NR(ER, mDM, sigma_nucleon):
    halo_model = hm.HaloModels()
    v_min = vmin_el(ER, mDM)
    v_max = hm.v_max(halo_model.v_esc)

    if v_min >= v_max:
        return 0

    def integrand(v):
        return (dsigmadER(ER, v, mDM, sigma_nucleon) * v * halo_model.velocity_dist_SHM(v))

    ret  = halo_model.rho_DM / mDM / mT() * quad(integrand,v_min, v_max)[0]

    return ret

#Differential rate per unit detector mass and recoil energy of elastic SI scattering
#ER:            recoil energy
#mDM:           mass of DM
#sigma_nucleon: DM-nucleon cross-section

def rateN_NR(Ne, mDM, sigma_nucleon):
    halo_model = hm.HaloModels()
    ER = LAr.Enr(Ne)
    v_min = vmin_el(ER, mDM)
    v_max = hm.v_max(halo_model.v_esc)
    rhoDM = halo_model.rho_DM

    if v_min >= v_max:
        return 0

    def integrand(v):
        return (dsigmadER(ER, v, mDM, sigma_nucleon) * v * halo_model.velocity_dist_SHM(v))

    ret  = LAr.dEnrdN(Ne) * rhoDM / mDM / mT() * quad(integrand,v_min, v_max)[0]

    return ret

#end
