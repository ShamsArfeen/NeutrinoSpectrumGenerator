#################################################
# ERWF.py file (Electron Radial Wave Fuunctions) Behrens Bühring description
#################################################
import numpy as np
from scipy.special import factorial2, binom

import IF
import Utilities
import Constants

#################################################
# Facteur d'expansion pour ERWF of electron : Hk, hk, Dk, dk
#################################################
def Hk(k, r, E, Z, R):
    sum = 0
    for mu in range(0, 1+1):
        #print('mu = ', mu)
        for nu in range(0, mu+1):
            #print('nu = ', nu)
            for rho in range(0, 2*nu+1):
                #print('rho = ', rho)
                if (rho== 0) :
                    I = 1
                elif (mu==1) and (nu==1) and (rho==1):
                    I = IF.I221(k, r, R)
                elif (mu==1) and (nu==1) and (rho==2):
                    I = IF.I222(k, r, R)
                else :
                    I = 0
                    print('problème avec I dans Hk', 'mu = ', mu, 'nu = ', nu, 'rho = ', rho)

                sum += factorial2(2*k-1, exact=True)/(factorial2(2*mu, exact=True)*factorial2(2*mu+2*k-1, exact=True)) * pow(-1, nu) * binom(mu, nu) * binom(2*nu, rho) * pow(r/R, 2*mu) * I * pow(R, 2*mu-2*nu) * pow(E*R, 2*nu-rho) * pow(Constants.ALPHA*Z, rho)

    Hk = sum
    return Hk

def hk(k, r, E, Z, R):
    sum = 0
    for mu in range(1, 1+1):
        #print('mu = ', mu)
        for nu in range(1, mu+1):
            #print('nu = ', nu)
            for rho in range(1, 2*nu):
                #print('rho = ', rho)
                I = IF.I211(k, r, R)
                sum = factorial2(2*k-1, exact=True)/(factorial2(2*mu, exact=True)*factorial2(2*mu+2*k-1, exact=True)) * pow(-1, nu) * binom(mu, nu) * binom(2*nu-1, rho) * pow(r/R, 2*mu) * I * pow(R, 2*mu-2*nu+1) * pow(E*R, 2*nu-1-rho) * pow(Constants.ALPHA*Z, rho)
    hk = sum
    return hk

def Dk(k, r, E, Z, R):
    sum = 0
    for mu in range(0, 1+1):
        #print('mu = ', mu)
        for nu in range(0, mu+1):
            #print('nu = ', nu)
            for rho in range(0, 2*nu+1+1):
                #print('rho = ', rho)
                if (rho==0) :
                    I = 1
                elif (mu==0) and (nu==0) and (rho==1):
                    I = IF.I111(k, r, R)
                elif (mu==1) and (nu==0) and (rho==1):
                    I = IF.I311(k, r, R)
                elif (mu==1) and (nu==1) and (rho==1):
                    I = IF.I331(k, r, R)
                elif (mu==1) and (nu==1) and (rho==2):
                    I = IF.I332(k, r, R)
                elif (mu==1) and (nu==1) and (rho==3):
                    I = IF.I333(k, r, R)
                else :
                    I = 0
                    print('problème avec I dans Dk', 'mu = ', mu, 'nu = ', nu, 'rho = ', rho)

                sum += factorial2(2*k-1, exact=True)/(factorial2(2*mu, exact=True)*factorial2(2*mu+2*k+1, exact=True)) * pow(-1, nu) * binom(mu, nu) * binom(2*nu+1, rho) * pow(r/R, 2*mu) * I * pow(R, 2*mu-2*nu) * pow(E*R, 2*nu+1-rho) * pow(Constants.ALPHA*Z, rho)

    Dk = sum
    return Dk

def dk(k, r, E, Z, R):
    sum = 0
    for mu in range(0, 1+1):
        #print('mu = ', mu)
        for nu in range(0, mu+1):
            #print('nu = ', nu)
            for rho in range(0, 2*nu+1):
                #print('rho = ', rho)
                if (rho==0) :
                    I = 1
                elif (mu==1) and (nu==1) and (rho==1):
                    I = IF.I321(k, r, R)
                elif (mu==1) and (nu==1) and (rho==2):
                    I = IF.I322(k, r, R)
                else :
                    I = 0
                    print('problème avec I dans dk', 'mu = ', mu, 'nu = ', nu, 'rho = ', rho)
                sum += factorial2(2*k-1, exact=True)/(factorial2(2*mu, exact=True)*factorial2(2*mu+2*k+1, exact=True)) * pow(-1, nu) * binom(mu, nu) * binom(2*nu, rho) * pow(r/R, 2*mu) * I * pow(R, 2*mu+1-2*nu) * pow(E*R, 2*nu-rho) * pow(Constants.ALPHA*Z, rho)
    dk = sum
    return dk

#################################################

#################################################
# ERWF of electron
#################################################
def fpk(k, r, E, Z, R):
    p = np.sqrt(pow(E, 2) - 1)
    fpk = Utilities.alphapk(k, E, Z, R) * pow(p*r, k-1)/factorial2(2*k-1) * (Hk(k, r, E) + hk(k, r, E))
    return fpk

def fmk(k, r, E, Z, R):
    p = np.sqrt(pow(E, 2) - 1)
    fmk = - Utilities.alphamk(k, E, Z, R) * pow(p*r, k-1)/factorial2(2*k-1) * (r/R) * (Dk(k, r, E) - dk(k, r, E))
    return fmk

def gpk(k, r, E, Z, R):
    p = np.sqrt(pow(E, 2) - 1)
    gpk = - Utilities.alphapk(k, E, Z, R) * pow(p*r, k-1)/factorial2(2*k-1) * (r/R) * (Dk(k, r, E) + dk(k, r, E))
    return gpk

def gmk(k, r, E, Z, R):
    p = np.sqrt(pow(E, 2) - 1)
    gmk = Utilities.alphamk(k, E, Z, R) * pow(p*r, k-1)/factorial2(2*k-1) * (Hk(k, r, E) - hk(k, r, E))
    return gmk

#################################################
