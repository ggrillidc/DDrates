#########################
#                       #
#     Migdal Effect     #
#                       #
#########################



################
# Explanations #
################

#This program contains the DM rate exploiting the Migdal Effect for Ar detectors,
#formulae from 1707.07258 and 1711.09906.

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
from scipy.integrate import quad, dblquad
import pandas as pd

import halo_models as hm
import DM_NR as DMNR
import DM_LAr as LAr


# Differential transition probabilities for Ar vs energy (eV), from the ancillary
# files of 1707.07258
df_migdal = pd.read_csv('../data/Migdal/migdal_transition_probability_Ar.csv')
# Relevant (n, l) electronic states
migdal_states = df_migdal.columns.values.tolist()
migdal_states.remove('E')

# Binding energies of the relevant Ar electronic states (in eV)
# From table II of 1707.07258
Enl_array = np.array([3.2e3,3.0e2, 2.4e2,27., 13.])#in eV
Enl_migdal = dict(zip(migdal_states,Enl_array))

#Minimum DM velocity to make a Migdal signal with energy E, with an elastic
#recoil ER and a DM mass mDM
def vmin_ME(E, ER, mDM):
    ynum = DMNR.mT() * ER + DMNR.mu_T(mDM) * E
    yden = DMNR.mu_T(mDM) * sqrt(2 * DMNR.mT() * ER)
    y = maximum(0, ynum / yden)
    return y


def rate_migdal(E, mDM, sigma_nucleon, nstate):

    halo_model = hm.HaloModels()
    ME_rate = 0
    for state, E_binding in Enl_migdal.items():
        E_binding *= nu.eV

        if nstate == 1:
            if state[0] not in ['1']:
                continue
        elif nstate == 2:
            if state[0] not in ['2']:
                continue
        elif nstate == 3:
            if state[0] not in ['3']:
                continue
        elif nstate == 12:
            if state[0] not in ['1','2']:
                continue

        prob = interp1d(df_migdal['E'].values * nu.eV,
                    df_migdal[state].values / nu.eV,
                    bounds_error=False,fill_value=0)

        def diff_rate(v, ER):
            Ee = E - E_binding

            if Ee < 0:
                return 0

            dsdER = DMNR.dsigmadER(ER, v, mDM, sigma_nucleon)

            SMHalo = halo_model.velocity_dist_SHM(v)

            #rescale the probability, see note added in 1707.07258
            qe2    = (nu.me * (2 * ER / DMNR.mT())**0.5 / (nu.eV / nu.c0))**2

            return (dsdER * v * SMHalo * qe2 / (2 * pi) * prob(Ee))

        r = dblquad(diff_rate, 0, DMNR.Emax(mDM, hm.v_max(halo_model.v_esc), DMNR.mT()),
            lambda ER: vmin_ME(E, ER, mDM), lambda _: hm.v_max(halo_model.v_esc))[0]

        ME_rate += r

    return halo_model.rho_DM / mDM / DMNR.mT() * ME_rate


def rate_migdalN(Ne, mDM, sigma_nucleon, nstate):

    E          = LAr.Eer(Ne)
    halo_model = hm.HaloModels()
    ME_rate = 0
    for state, E_binding in Enl_migdal.items():
        E_binding *= nu.eV

        if nstate == 1:
            if state[0] not in ['1']:
                continue
        elif nstate == 2:
            if state[0] not in ['2']:
                continue
        elif nstate == 3:
            if state[0] not in ['3']:
                continue
        elif nstate == 12:
            if state[0] not in ['1','2']:
                continue

        prob = interp1d(df_migdal['E'].values * nu.eV,
                        df_migdal[state].values / nu.eV,
                        bounds_error=False,fill_value=0)

        def diff_rate(v, ER):
            Ee = E - E_binding

            if Ee < 0:
                return 0

            dsdER = DMNR.dsigmadER(ER, v, mDM, sigma_nucleon)

            SMHalo = halo_model.velocity_dist_SHM(v)
            #rescale the probability, see note added in 1707.07258
            qe2    = (nu.me * (2 * ER / DMNR.mT())**0.5 / (nu.eV / nu.c0))**2

            return (dsdER * v * SMHalo * qe2 / (2 * pi) * prob(Ee))

        r = dblquad(diff_rate, 0, DMNR.Emax(mDM, hm.v_max(halo_model.v_esc), DMNR.mT()),
            lambda ER: vmin_ME(E, ER, mDM), lambda _: hm.v_max(halo_model.v_esc))[0]

        ME_rate += r

    rateE = halo_model.rho_DM / mDM / DMNR.mT() * ME_rate
    rateN = LAr.dEerdN(Ne)*rateE

    return rateN

#end
