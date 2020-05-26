#####################################
#                                   #
#     DarkSide useful functions     #
#                                   #
#####################################



################
# Explanations #
################

#This program contains some useful functions for analyses of LAr experiments,
#from 1802.06998 and 1802.06994


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



dataER = np.genfromtxt("../data/DarkSide50/NeVSEer2.txt", delimiter=' ')
x1 = dataER[:,0]
y1 = dataER[:,1]
Nel1 = interp1d(x1, y1, fill_value='extrapolate')
Eer1 = interp1d(y1, x1, fill_value='extrapolate')

#Calibration curve used to convert electron recoil spectra to ionization spectra,
#from Fig. 2 of 1802.06998.
#Eer:        electron recoil energy in keV
def NeEr(Eer):
    if Eer < 0.17177770796823463:
        ret = 0.103611 + 45.9686 * Eer
    elif Eer <= 3.:
        ret = Nel1(Eer)
    else:
        ret = 7.150456384177961 * Eer + 27.849112071870415
    return ret


#Calibration curve used to convert electron recoil spectra to ionization spectra,
#from Fig. 2 of 1802.06998.
#Ne:        number of electrons

def Eer(Ne):
    if Ne <= 8:
        ret = -0.021754 * (0.103611 - Ne) * nu.keV
    elif Ne <= 49:
        ret = Eer1(Ne)*nu.keV
    else:
        ret = -0.139851 * (27.8491 - Ne) * nu.keV
    return ret


#Derivative of the calibration curve used to convert electron recoil spectra to ionization spectra,
#from Fig. 2 of 1802.06998. Fit obtained with mathematica.
#This is the Jacobian of the transformation dR/dE to dR/dNe
#Ne:        number of electrons

def dEerdN(Ne):
    a2 = 0.025044254616607853
    a3 = 0.000023491507516026288
    a4 = 0.000014616419760867523
    if Ne <= 8:
        ret = 0.021754 * nu.keV
    elif Ne <= 49:
        ret = (a2 + 2 * a3 * Ne + 3 * a4 * Ne**2) * nu.keV
    else:
        ret = 0.139851 * nu.keV
    return ret


#Calibration curve used to convert nuclear recoil spectra to ionization spectra,
#obtained from Fig. 6 of 1802.06994.
dataNR = np.genfromtxt("../data/DarkSide50/quenchingDS.txt", delimiter=' ')

x = dataNR[:,1]
y = dataNR[:,0]
Enr1 = interp1d(x, y, fill_value='extrapolate')
Nnr1 = interp1d(y, x, fill_value='extrapolate')

def Enr(Ne):
    return Enr1(Ne)* nu.keV

def Nnr(Enr):
    return Nnr1(Enr)

#Derivative of the calibration curve used to convert nuclear recoil spectra to ionization spectra,
#from Fig. 6 of 1802.06994.
#This is the Jacobian of the transformation dR/dE to dR/dNe
#Ne:        number of electrons

def dEnrdN(Ne):
    eps = 0.1
    df = (Enr(Ne+eps) - Enr(Ne-eps)) / (2 * eps)
    return df

#end
