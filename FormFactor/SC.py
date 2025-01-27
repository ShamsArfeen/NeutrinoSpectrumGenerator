#################################################
# SC.py file (Shape Factor)
#################################################
import numpy as np
import matplotlib.pyplot as plt

import FF
import Constants
import MK_mK
import Utilities


E0 = Constants.E0
Z = Constants.Z

def shape(E, chemin):
    p = np.sqrt(pow(E, 2) - 1)
    q = E0 - p
    R = 1.2*pow(210, 1/3) / 385
    fv0101_0 = 0.0465 # FF.FVNKLs(0, 1, 0, 1, 1, 0, 0, 0, chemin)

    fv0110_1 = FF.FVNKLs(0, 1, 1, 0, 1, 1, 1, 1, chemin)
    fv0110_0 = FF.FVNKLs(0, 1, 1, 0, 1, 0, 0, 0, chemin)
    fa0111_0 = FF.FANKLs(0, 1, 1, 1, 1, 0, 0, 0, chemin)
    fa0111_1 = FF.FANKLs(0, 1, 1, 1, 1, 1, 1, 1, chemin)
    
    print('fv0101_0 : ', fv0101_0)
    print('fv0110_1 : ', fv0110_1)
    print('fv0110_0 : ', fv0110_0)
    print('fa0111_0 : ', fa0111_0)
    print('fa0111_1 : ', fa0111_1)
    #print('mu ke : ', Utilities.mu_ke(1, E, chemin))
    #print('lambda ke=1 : ', Utilities.lambda_ke(1, E, chemin))
    #print('lambda ke=2 : ', Utilities.lambda_ke(2, E, chemin))
    #print('gamma ke : ', Utilities.gamma_ke(1, chemin))
    
    M111 = - fv0101_0 - pow(1/3, 1.5)*Constants.ALPHA*Z*fv0110_1  - pow(1/3, 1.5)*E0*R*fv0110_0 - 1/3*Constants.ALPHA*Z*np.sqrt(2/3)*fa0111_1 - 1/3*(E - q)*R*np.sqrt(2/3)*fa0111_0
    m111 = -1/3*R*(np.sqrt(1/3)*fv0110_0 + np.sqrt(2/3)*fa0111_0)
    
    M121 = 1/3*p*R*(np.sqrt(2/3)*fv0110_0 - np.sqrt(1/3)*fa0111_0 )
    M112 = 1/3*q*R*(np.sqrt(2/3)*fv0110_0 + np.sqrt(1/3)*fa0111_0)
    
    m121 = MK_mK.Kn(1, 2, 1)*p*R*R/5*(-fv0110_0 - np.sqrt(2)*fa0111_0)
    m112 = MK_mK.Kn(1,1,2)*q*R*R/3*(-fv0110_0 - np.sqrt(2)*fa0111_0)
    
    #print('M111 : ', M111)
    s = Utilities.lambda_ke(1, E, chemin)*(M111**2 + m111**2 - 2*Utilities.mu_ke(1, E, chemin)*Utilities.gamma_ke(1, chemin)/E * M111*m111) \
    + Utilities.lambda_ke(1, E, chemin)*(M111**2-2*Utilities.mu_ke(1, E, chemin)*Utilities.gamma_ke(1, chemin)/E * M111*m111) \
    + Utilities.lambda_ke(2, E, chemin)*(M121**2-2*Utilities.mu_ke(2, E, chemin)*Utilities.gamma_ke(2, chemin)/2/E*M121*m121) \
    + Utilities.lambda_ke(1, E, chemin)*(M112**2-2*Utilities.mu_ke(1, E, chemin)*Utilities.gamma_ke(1, chemin)/E*M112*m112)
    
    return s

def shape_B(E, chemin):
    p = np.sqrt(pow(E, 2) - 1)
    q = E0 - p
    R = 1.2*pow(210, 1/3) / 385
    fv0101_0 = 0.0465 # 0.0515 # 
    fv0110_1 = -0.1438 # -0.1694 # 
    fv0110_0 = -0.1108 # -0.1413 # 
    fa0111_0 = -0.1776 # -0.1993 # 
    fa0111_1 = -0.1580 # -0.1722 # 

    
    M111 = - fv0101_0 - pow(1/3, 1.5)*Constants.ALPHA*Z*fv0110_1  - pow(1/3, 1.5)*E0*R*fv0110_0 - 1/3*Constants.ALPHA*Z*np.sqrt(2/3)*fa0111_1 - 1/3*(E - q)*R*np.sqrt(2/3)*fa0111_0
    m111 = -1/3*R*(np.sqrt(1/3)*fv0110_0 + np.sqrt(2/3)*fa0111_0)
    
    M121 = 1/3*p*R*(np.sqrt(2/3)*fv0110_0 - np.sqrt(1/3)*fa0111_0 )
    M112 = 1/3*q*R*(np.sqrt(2/3)*fv0110_0 + np.sqrt(1/3)*fa0111_0)
    
    m121 = MK_mK.Kn(1, 2, 1)*p*R*R/5*(-fv0110_0 - np.sqrt(2)*fa0111_0)
    m112 = MK_mK.Kn(1,1,2)*q*R*R/3*(-fv0110_0 - np.sqrt(2)*fa0111_0)
    
    #print('M111 : ', M111)
    s = Utilities.lambda_ke(1, E, chemin)*(M111**2 + m111**2 - 2*Utilities.mu_ke(1, E, chemin)*Utilities.gamma_ke(1, chemin)/E * M111*m111) + Utilities.lambda_ke(1, E, chemin)*(M111**2-2*Utilities.mu_ke(1, E, chemin)*Utilities.gamma_ke(1, chemin)/E * M111*m111) + Utilities.lambda_ke(2, E, chemin)*(M121**2 - 2*Utilities.mu_ke(2, E, chemin)*Utilities.gamma_ke(2, chemin)/2/E*M121*m121)+Utilities.lambda_ke(1, E, chemin)*(M112**2-2*Utilities.mu_ke(1, E, chemin)*Utilities.gamma_ke(1, chemin)/E*M112*m112)
    
    return s

            
chemin_bi = "./test-Bi215.txt"

NUM_BINS = 4500
E_MAX = 4.5 # MeV
n = int(NUM_BINS * (E0 - 1) / (E_MAX/0.51099895069))

energy = np.linspace(1, E0, n)

y = shape(energy, chemin_bi)
print(y.shape)
np.save('E_Bi-215.npy', energy)
np.save('Bi-215.npy', y)

plt.figure('')
plt.plot(energy, shape(energy, chemin_bi), label='Computed C(Z,W)')
# plt.plot(energy, shape_B(energy, chemin_bi), label='Behrens, 1970')
plt.xlabel('W')
plt.ylabel(r'$C(E)$')
plt.title(f"Bi-215, GS (9/2-) to GS (9/2+) decay shape factor, W0 = {str(E0)}")
plt.legend()
plt.grid()
plt.savefig('Shape Factor for 215Bi')
plt.show()
