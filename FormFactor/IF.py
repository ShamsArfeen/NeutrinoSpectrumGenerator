#################################################
# IF.py file (I Function)
#################################################
import numpy as np
from scipy import integrate

import Utilities

#################################################
# I(k, m, n, rho, r) Functions
#################################################
def I111(k, r, R):
    if (r >= 0) and (r <= R) :
        i111 = (3/2) - (2*k+1)/(2*(2*k+3)) * pow(r/R, 2)
    elif r > R : 
        i111 = (2*k+1)/(2*k) * (R/r) - (3/(2*k*(2*k+3))) * pow(R/r, 2*k+1)
    return i111
   
def I222(k, r, R):
    if (r >= 0) and (r <= R) :
        i222 = (9/4) - (3*(k+1))/(2*(2*k+3)) * pow(r/R, 2) + (2*k+1)/(12*(2*k+3))* pow(r/R, 4)
    elif r > R : 
        i222 = (2*k+1)/(k) * pow(R/r, 2)*np.log(r/R) + ( 19/12 + 7/(12*(2*k+3)) - 3/(2*pow(k, 2)*(2*k+3)) )*pow(R/r, 2) + 3/(2*pow(k, 2)*(2*k+3)) * pow(R/r, 2*k+2)
    return i222

def I221(k, r, R):
    if (r >= 0) and (r <= R) :
        i221 = (3/2) - (k+1)/(2*(2*k+3)) * pow(r/R, 2)
    elif r > R : 
        i221 = (4*k+1)/(2*k) * (R/r) - (3*k)/(2*(2*k-1))*pow(R/r, 2) + 3/(2*k*(2*k+3)*(2*k-1)) * pow(R/r, 2*k+1)
    return i221

def I211(k, r, R):
    if (r >= 0) and (r <= R) :
        i211 = - 1/(2*(2*k+3)) * pow(r/R, 2)
    elif r > R : 
        i211 = - (1/k) * (R/r) + 3/(2*(2*k-1))*pow(R/r, 2) - 3/(k*(2*k+3)*(2*k-1)) * pow(R/r, 2*k+1)
    return i211

def I333(k, r, R):
    if (r >= 0) and (r <= R) :
        i333 = 27/8 - (9*(4*k+5))/(8*(2*k+5))* pow(r/R, 2) + (8*k+7)/(8*(2*k+7))* pow(r/R, 4) - (2*k+1)/(24*(2*k+9))* pow(r/R, 6)
    elif r > R : 
        i333 = (2*k+1)*(2*k+3)/(2*pow(k, 2)) * pow(R/r, 3)*np.log(r/R) + (19/12 + 5/(3*k) - 2/pow(k, 2) - 3/(2*pow(k, 3)))*pow(R/r, 3) + 3/(2*pow(k, 2))*pow(R/r, 2*k+3)*np.log(r/R) + (45/(8*(2*k+5)) - 21/(8*(2*k+7)) + 1/(3*(2*k+9)) - 5/(3*k) + 2/pow(k, 2) + 3/(2*pow(k, 2)))*pow(R/r, 2*k+3)
    return i333

def I332(k, r, R):
    if (r >= 0) and (r <= R) :
        i332 = 9/4 - (4*k+5)/(2*(2*k+5))* pow(r/R, 2) + (8*k+7)/(36*(2*k+7))* pow(r/R, 4)
    elif r > R : 
        i332 = (2*k+3)/(3*k)* pow(R/r, 2)*np.log(r/R) + (67/36 + 23/(12*(2*k+1)) + 1/k - 1/(2*pow(k, 2)))* pow(R/r, 2) - (2*k+3)/(2*(2*k-1))* pow(R/r, 3) + (4*k-1)/(2*pow(k, 2)*(2*k-1))* pow(R/r, 2*k+2) + (-23/(12*(2*k+1)) + 5/(2*(2*k+5)) - 7/(12*(2*k+7)))* pow(R/r, 2*k+3)
    return i332

def I322(k, r, R):
    if (r >= 0) and (r <= R) :
        i322 = 9/4 - (3*(2*k+3))/(4*(2*k+5))*pow(r/R, 2) + (k+2)/(6*(2*k+7))*pow(r/R, 4)
    elif r > R : 
        i322 = (2*k+3)/(k)* pow(R/r, 2)*np.log(r/R) + (19/12 + 69/(12*(2*k+1)) - 3/k - 3/(2*pow(k, 2)))* pow(R/r, 2) + (3*(2*k+3))/(4*k*(2*k-1))* pow(R/r, 3) - 3/(2*pow(k, 2)*(2*k-1))* pow(R/r, 2*k+2) + (3/(2*(2*k+5)) - 1/(4*(2*k+7)) - 23/(4*(2*k+1)) + 9/(4*k))* pow(R/r, 2*k+3)
    return i322

def I331(k, r, R):
    if (r >= 0) and (r <= R) :
        i331 = 3/2 - (4*k+5)/(6*(2*k+5)) * pow(r/R, 2)
    elif r > R : 
        i331 = (5*k+1)*(2*k+3)/(6*k*(k+1))* (R/r) - (k*(2*k+3))/((2*k-1)*(2*k+1))* pow(R/r, 2) + 1/(2*k*(2*k-1))* pow(R/r, 2*k+1) + (7/6 + 5/(6*(2*k+5)) - (5*k+1)*(2*k+3)/(6*k*(k+1)) + (k*(2*k+3))/((2*k-1)*(2*k+1)) - 1/(2*k*(2*k-1)))* pow(R/r, 2*k+3)
    return i331

def fct1(x, Z, R):
    return x * Utilities.Ur(x, Z, R)

def fct2(x, k, Z, R):
    return pow(x, 2 * k + 2) * Utilities.Ur(x, Z, R)

def fct3(x, k, Z, R):
    return pow(x, 2 * k) * Utilities.Ur(x, Z, R)

def I311(k, r, Z, R):
 
    resultat1 = integrate.quad(fct1, 0, r, args=(Z, R))
    resultat2 = integrate.quad(fct2, 0, r, args=(k, Z, R))
    resultat3 = integrate.quad(fct3, 0, r, args=(k, Z, R))
        
    i311 = ( (4*(k+1)*(2*k+3))/(2*k+1) * pow(r, -2*k-3) * resultat2[0] - ((k+1)*(2*k+3))/(2*k+1) * pow(r, -2*k-1) * resultat3[0] + (4*(2*k+3))/((2*k+1)*(2*k-1)) * pow(r, -2) * resultat1[0] )
    return i311

def I321(k, r, Z, R):
    
    resultat1 = integrate.quad(fct1, 0, r, args=(Z, R))
    resultat2 = integrate.quad(fct2, 0, r, args=(k, Z, R))
    i321 = ((2*(2*k+3))/(2*k+1) * pow(r, -2) * resultat1[0] - (2*(2*k+3))/(2*k+1) * pow(r, -2*k-3) * resultat2[0])
    return i321  

def I(k, r, m, n, rho, R):
    I = 0
    if rho==0 :
        I = 1
    elif (m==1) and (n==1) and (rho==1):
        I = I111(k, r, R)
    elif (m==2) and (n==1) and (rho==1):
        I = I211(k, r, R)
    elif (m==2) and (n==2) and (rho==1):
        I = I221(k, r, R)
    elif (m==2) and (n==2) and (rho==2):
        I = I222(k, r, R)
    elif (m==3) and (n==3) and (rho==1):
        I = I331(k, r, R)
    elif (m==3) and (n==2) and (rho==2):
        I = I322(k, r, R)
    elif (m==3) and (n==3) and (rho==2):
        I = I332(k, r, R)
    elif (m==3) and (n==3) and (rho==3):
        I = I333(k, r, R)
    return I
#################################################
