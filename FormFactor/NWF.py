#################################################
# NWF.py file (Nuclear Wave Functions)
#################################################
import numpy as np
import math
from scipy.special import sph_harm, binom, eval_genlaguerre
from sympy.physics.quantum.cg import CG
from sympy.physics.wigner import wigner_9j
from sympy import S

import Constants


#################################################
# Nuclear wave function with normalization calculation
#################################################
def factorial(n):
    return math.factorial(n)

def double_factorial(n):
    if n <= 0:
        return 1
    return np.prod(range(n, 0, -2))

def radial_ho(n, l, nu, r):
    k = n - 1
    
    # Facteur de normalisation
    N = np.sqrt(
        np.sqrt(2 * nu**3 / np.pi) * 2**(k + 2 * l + 3) *
        factorial(k) * nu**l /
        double_factorial(2 * k + 2 * l + 1)
    )

    # Polynôme de Laguerre généralisé
    gl = eval_genlaguerre(k, l + 0.5, 2*nu * r**2)
    wf = N * r**(l) * np.exp(-nu * r**2) * gl

    return wf

def get_occupation_numbers(N):
    occ_numbers = [1, 0, 1, min(N, 2)]

    for i in range(4, len(Constants.sm_occupation), 4):
        if N > Constants.sm_occupation[i - 1]:
            occ_numbers.append(Constants.sm_occupation[i])
            occ_numbers.append(Constants.sm_occupation[i + 1])
            occ_numbers.append(Constants.sm_occupation[i + 2])
            occ_numbers.append(min(N - Constants.sm_occupation[i - 1], 
                                    Constants.sm_occupation[i + 3] - Constants.sm_occupation[i - 1]))
    return occ_numbers
    
def calc_nu(rms, Z):
    occ_numbers = get_occupation_numbers(Z)

    s = 0.0
    for i in range(0, len(occ_numbers), 4):
        s += (4 * occ_numbers[i] + 2 * occ_numbers[i + 1] - 1) * occ_numbers[i + 3]

    return s / Z / 2.0 / rms / rms

def NWF(n, l, nu, r): # l = k (kappa)
    if l == -1 :
        l = 0
        
    dr = r[1]-r[0]
    LNWF = radial_ho(n, l, nu, r)
    SNWF = []
    for i in range(len(r)-1):
        SNWF.append(1/(2*Constants.UMASSC2/Constants.me) * ((LNWF[i+1]-LNWF[i])/dr + (l+1)/r[i] * LNWF[i]))     
    SNWF.append( SNWF[-1] )
    
    Lunity = 0
    Sunity = 0
    for i in range(len(r)-1):
        Lunity += (pow(LNWF[i+1]*r[i+1], 2) + pow(LNWF[i]*r[i], 2))/2 * dr
        Sunity += (pow(SNWF[i+1]*r[i+1], 2) + pow(SNWF[i]*r[i], 2))/2 * dr
    unity = Lunity + Sunity
    #print('lunity : ', Lunity)
    #print('sunity : ', Sunity)    
    return LNWF, SNWF#/np.sqrt(unity)

#################################################

#################################################
# Spherical Harmonique (scalare and vector)
#################################################
def Y(M, L, teta, psi):
    
    #teta = np.linspace(0, np.pi, 200)
    #psi = np.linspace(0, 2*np.pi, 200)
    #teta, psi = np.meshgrid(teta, psi)
    #xyz = np.array([np.sin(teta)*np.sin(psi), np.sin(teta)*np.cos(psi), np.cos(teta)])
    
    Y = sph_harm(M, L, teta, psi)
    if M < 0 :
        Y = np.sqrt(2) * pow(-1, M) * np.imag(Y)
    elif M > 0:
        Y = np.sqrt(2) * pow(-1, M) * np.real(Y)
    else :
        Y = sph_harm(0, L, teta, psi)

    #Yx, Yy, Yz = xyz * abs(Y)
    return Y    
    
def TMLL0(M, L, teta, psi):
    return pow(1j, L) * Y(M, L, teta, psi)

def TMKL1(M, K, L, teta, psi, j1, m1, j2, j3):
    ep1 = - 1/np.sqrt(2) * np.array([1, 1j, 0])
    em1 = 1/np.sqrt(2) * np.array([1, -1j, 0])
    e0 = np.array([0, 0, 1])
    j2 = 2

    T = (CG(S(j1)/2, S(m1)/2, S(j2)/2, S(-2)/2, S(j3)/2, S(2)/2).doit()*ep1*Y(M, L, teta, psi) + CG(S(j1)/2, S(m1)/2, S(j2)/2, S(0)/2, S(j3)/2, S(0)/2).doit()*e0*Y(M, L, teta, psi) + CG(S(j1)/2, S(m1)/2, S(j2)/2, S(2)/2, S(j3)/2, S(-2)/2).doit()*em1*Y(M, L, teta, psi))
    return T * pow(1j, L) * pow(-1, L-K+1)
#################################################
