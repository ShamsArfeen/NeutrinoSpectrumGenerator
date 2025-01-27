#################################################
# Utilities.py file
#################################################
import numpy as np
from scipy.special import factorial2, gamma
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_9j
from sympy import S
import matplotlib.pyplot as plt

import Constants
#################################################
# Read OBTDs file
#################################################
def find_nucleus(chemin):
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        lines = fichier_txt.readlines()
        beta_type = str(lines[0].split()[1])
        A = int(lines[0].split()[2])
        Z = int(lines[0].split()[3])
        Ji = float(lines[0].split()[4])
        nucleus_name = str(lines[0].split()[5])
        f.close()
    return beta_type, A, Z, Ji, nucleus_name
#################################################

#################################################
# Radius
#################################################
def calc_R(A) :
    return (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
#################################################

#################################################
# End point energy
#################################################
def end_point(mi, mf, betatype) :
    delta = mi - mf
    if betatype == 'beta-' :
        end_point = (delta + Constants.me)/Constants.me
    elif betatype == 'beta+' :
        end_point = (delta - 2*Constants.me)/Constants.me
    Ecp_max = ( (((delta-Constants.me)/Constants.me)**2) + 2 * ((delta-Constants.me)/Constants.me))/(2*(mf/Constants.me))
    return end_point - Ecp_max
#################################################

#################################################
# Fermi Function
#################################################
def fermi_function(W, Z, R):
    """Traditional Fermi Function
    :param Z: Proton number of the final nuclear state
    :param W: Electron energy in units of me c^2
    :param R: Nuclear radius in units of the electron Compton wavelength
    """
    f = 1

    if Z == 0:
        return f

    g = np.sqrt(1-(Constants.ALPHA*Z)**2)
    p = np.sqrt(W**2-1)
    y = Constants.ALPHA*Z*W/p

    #We use the traditional Fermi function, i.e. a prefactor 4 instead of 2(1+gamma)
    #This is consistent with the L0 description below

    f = (4
            *np.power(2*p*R, 2*(g-1))
            *np.exp(np.pi*y)
            *(np.abs(gamma(g+1.j*y))/(gamma(1+2*g)))**2)

    return f

def GenFermiFunction(k, E, Z, R):
    p = np.sqrt(pow(E, 2) - 1)
    g = np.sqrt(k-(Constants.ALPHA*Z)**2)
    y = Constants.ALPHA*Z*E/p

    fermi = pow(k*factorial2(2*k-1), 2)*pow(4, k)*pow(2*p*R, 2*(g-k))*np.exp(np.pi*y)*pow(abs(gamma(g+1j*y))/gamma(1+2*g), 2)
    return fermi

#################################################

#################################################
# Potential U(r) and V(r)
#################################################
def Ur(r, Z, R):
    if r <= R :
        return 3/2 - 1/2 * pow(r/R, 2)
    elif r > R :
        return R/r

def V(r, Z, R):
    v = 1.44*Z/R * Ur(r, Z, R)
    return v
#################################################

#################################################
# Beta Spectrum and integration
#################################################
def beta_spectrum(E, Z, R, mi, mf, betatype):
    E0 = end_point(mi, mf, betatype)
    p = np.sqrt(pow(E, 2) - 1)
    spectrum = pow(p, 2)*pow(E0 - E, 2)*fermi_function(E, Z, R)
    return spectrum

def neutrino_spectrum(E, Z, R, mi, mf, betatype):
    E0 = end_point(mi, mf, betatype)
    p = np.sqrt(pow(E, 2) - 1)
    pnu =  E0 - p
    Enu = pnu
    
    spectrum = pow(pnu, 2) * pow(E0 - Enu, 2)#*fermi_function(E, Z, R)
    return spectrum
    
"""mi = 209.984120237*Constants.UMASSC2
mf = 209.982873686*Constants.UMASSC2
energy = np.linspace(1, end_point(mi, mf, 'beta-'), 1000)
plt.figure('')
plt.plot(energy, beta_spectrum(energy, 83, 1.27*pow(210, 1/3), mi, mf, 'beta-'), label='beta spectrum')
plt.plot(energy, neutrino_spectrum(energy, 83, 1.27*pow(210, 1/3), mi, mf, 'beta-'), label='neutrino spectrum')
plt.legend()
plt.grid()
plt.show()"""
    
def integrate_betaspectrum(Z, R):
    result = integrate.quad(betaspectrum, 1, delta, args=(Z, R))
    return result[0]
#################################################

#################################################
# Coulomb Function
#################################################
def alpha_k(k, E, Z, R): # eq to beta in Banbynek paper (not sure)
    p = np.sqrt(pow(E, 2) - 1)
    g = np.sqrt(abs(k)-(Constants.ALPHA*Z)**2)
    alphak = np.sqrt(GenFermiFunction(abs(k), E, Z, R))*pow(abs(k), -1)*p*np.sqrt(2*E)*np.sqrt((abs(k)+g)*(abs(k)*E-np.sign(k)*g))
    return alphak

def gamma_ke(ke, chemin):
    beta_type, A, Z_mother, _, _ = find_nucleus(chemin)
    R = 1.2*pow(A, 1/3)
    if beta_type == 'beta+' :
        Z = Z_mother - 1
    else :
        Z = Z_mother + 1
    gke = np.sqrt(pow(ke, 2) - pow(Constants.ALPHA * Z, 2))
    return gke

def lambda_ke(ke, E, chemin):
    beta_type, A, Z_mother, _, _ = find_nucleus(chemin)
    R = 1.2*pow(A, 1/3)
    if beta_type == 'beta+' :
        Z = Z_mother - 1
    else :
        Z = Z_mother + 1
    #lke = (pow(alpha_k(-ke, E, Z, R), 2) + pow(alpha_k(ke, E, Z, R), 2))/(pow(alpha_k(-1, E, Z, R), 2) + pow(alpha_k(1, E, Z, R), 2))
    gk = gamma_ke(ke, chemin)
    g1 = gamma_ke(1, chemin)
    lke = GenFermiFunction(ke, E, Z, R)/GenFermiFunction(1, E, Z, R)*(ke+gk)/(ke*(1+g1))
    return lke

def mu_ke(ke, E, chemin):
    beta_type, A, Z_mother, _, _ = find_nucleus(chemin)
    R = 1.2*pow(A, 1/3)
    if beta_type == 'beta+' :
        Z = Z_mother - 1
    else :
        Z = Z_mother + 1

    mke = ke * E/gamma_ke(ke, chemin) * (pow(alpha_k(-ke, E, Z, R), 2) - pow(alpha_k(ke, E, Z, R), 2))/(pow(alpha_k(-ke, E, Z, R), 2) + pow(alpha_k(ke, E, Z, R), 2))
    return mke  
#################################################

#################################################
# Factor GKLs
#################################################
def GKLs(K, L, s, kf, ki):
    lnf = kf
    lni = ki
    jf = lnf - 1/2
    ji = lni - 1/2
    if s == 0 :
        L = K
        
    if lnf < 0 :
        lnf = abs(lnf) - 1
        jf = lnf + 1/2
    if lni < 0 : 
        lni = abs(lni) - 1
        ji = lni + 1/2

    CGcoeff = CG(S(lnf*2)/2, S(0)/2, S(lni*2)/2, S(0)/2, S(L*2)/2, S(0)/2).doit()
    WC = wigner_9j(K, s, L, jf, 1/2, lnf, ji, 1/2, lni, prec=64).doit()
    GKLs = np.sqrt((2*s+1) * (2*K+1) * (2*lnf+1) * (2*lni+1) * (2*jf+1) * (2*ji+1)) * pow(1j, lni+lnf+L) * pow(-1, ji-jf) * CGcoeff * WC 
    return float(GKLs)
#################################################
