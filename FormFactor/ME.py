#################################################
# ME.py file (Matrix Elements)
#################################################
import numpy as np

import NWF
import IF
import Utilities
import Constants


#################################################
# Matrix Elements
#################################################
def MVNKK0(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R):    
    lnf = int(kf)
    lni = int(ki)

    if lnf < 0 :
        lnf = abs(lnf) - 1
    if lni < 0 : 
        lni = abs(lni) - 1
    
    #b = 0.836 * R * pow(210, -1/6) 
    #nui = nuf = 1/(2*pow(b, 2))
    gf, ff = NWF.NWF(nf, lnf, nuf, r)
    gi, fi = NWF.NWF(ni, lni, nui, r)

    II = []
    for i in range(len(r)):
        II.append(IF.I(k, r[i], m, n, rho, R))

    fct1 = gf*pow(r/R, K+2*N)*II*gi*pow(r, 2)
    fct2 = ff*pow(r/R, K+2*N)*II*fi*pow(r, 2)
  
    int1, int2 = 0, 0
    dr = r[1] - r[0]
    for i in range(len(r)-1):
        int1 += (fct1[i]+fct1[i+1])/2 * dr
        int2 += (fct2[i]+fct2[i+1])/2 * dr
    MVNKK0 = np.sqrt(2)/np.sqrt(2*Ji+1) * (Utilities.GKLs(K, K, 0, kf, ki) * int1 + np.sign(kf)*np.sign(ki) * Utilities.GKLs(K, K, 0, -kf, -ki) * int2)
    
    return MVNKK0 * Constants.gv, fct1, fct2

def MANKK0(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R):    
    lnf = int(kf)
    lni = int(ki)

    if lnf < 0 :
        lnf = abs(lnf) - 1
    if lni < 0 : 
        lni = abs(lni) - 1
        
    gf, ff = NWF.NWF(nf, lnf, nuf, r)*pow(r, 2)
    gi, fi = NWF.NWF(ni, lni, nui, r)*pow(r, 2)

    II = []
    for i in range(len(r)):
        II.append(IF.I(k, r[i], m, n, rho, R))

    fct1 = gf*pow(r/R, K+2*N)*II*fi
    fct2 = ff*pow(r/R, K+2*N)*II*gi

    int1, int2 = 0, 0
    dr = r[1] - r[0]
    for i in range(len(r)-1):
        int1 += (fct1[i]+fct1[i+1])/2 * dr
        int2 += (fct2[i]+fct2[i+1])/2 * dr

    MANKK0 = np.sqrt(2)/np.sqrt(2*Ji+1) * (np.sign(ki)* Utilities.GKLs(K, L, 0, kf, -ki) * int1 + np.sign(kf) * Utilities.GKLs(K, L, 0, -kf, ki) * int2)
    return MANKK0, fct1, fct2

def MVNKL1(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R):
    lnf = int(kf)
    lni = int(ki)

    if lnf < 0 :
        lnf = abs(lnf) - 1
    if lni < 0 : 
        lni = abs(lni) - 1
        
    gf, ff = NWF.NWF(nf, lnf, nuf, r)
    gi, fi = NWF.NWF(ni, lni, nui, r)

    II = []
    for i in range(len(r)):
        II.append(IF.I(k, r[i], m, n, rho, R))
        
    fct1 = gf*pow(r/R, L+2*N)*II*fi*pow(r, 2)
    fct2 = ff*pow(r/R, L+2*N)*II*gi*pow(r, 2)

    int1, int2 = 0, 0
    dr = r[1] - r[0]
    for i in range(len(r)-1):
        int1 += (fct1[i]+fct1[i+1])/2 * dr
        int2 += (fct2[i]+fct2[i+1])/2 * dr
        
    MVNKL1 = np.sqrt(2)/np.sqrt(2*Ji+1) * (np.sign(ki)*Utilities.GKLs(K, L, 1, kf, -ki) * int1 + np.sign(kf) * Utilities.GKLs(K, L, 1, -kf, ki) * int2)
    return MVNKL1, fct1, fct2

def MANKL1(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R):  
    lnf = int(kf)
    lni = int(ki)

    if lnf < 0 :
        lnf = abs(lnf) - 1
    if lni < 0 : 
        lni = abs(lni) - 1
    
    #b = 0.836 * R * pow(210, -1/6) 
    #nui = nuf = 1/(2*pow(b, 2))    
    gf, ff = NWF.NWF(nf, lnf, nuf, r)
    gi, fi = NWF.NWF(ni, lni, nui, r)

    II = []
    for i in range(len(r)):
        II.append(IF.I(k, r[i], m, n, rho, R))
    
    fct1 = gf*pow(r/R, L+2*N)*II*gi*pow(r, 2)
    fct2 = ff*pow(r/R, L+2*N)*II*fi*pow(r, 2)

    int1, int2 = 0, 0
    dr = r[1] - r[0]
    for i in range(len(r)-1):
        int1 += (fct1[i]+fct1[i+1])/2 * dr
        int2 += (fct2[i]+fct2[i+1])/2 * dr
    MANKL1 = np.sqrt(2)/np.sqrt(2*Ji+1) * (Utilities.GKLs(K, L, 1, kf, ki) * int1 + np.sign(kf)*np.sign(ki) * Utilities.GKLs(K, L, 1, -kf, -ki) * int2)
    
    return MANKL1, fct1, fct2
    
def MVNKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R):
    if s==0 :
        MVNKLs = MVNKK0(N, K, L, 0, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R)[0]
    elif s==1 :
        MVNKLs = MVNKL1(N, K, L, 1, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R)[0]
    return float(MVNKLs)

def MANKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R):
    if s==0 :
        MANKLs = MANKK0(N, K, L, 0, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R)[0]
    elif s==1 :
        MANKLs = MANKL1(N, K, L, 1, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R)[0]
    return float(MANKLs)
