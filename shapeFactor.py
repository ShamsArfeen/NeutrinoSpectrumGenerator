
import sys

sys.path.append("./FormFactor")

from FormFactor.FF import FANKLs, FVNKLs
import FormFactor.Utilities as Utilities
from Resolution import *

import numpy as np

Me = 0.51099895069
ALPHA = 0.0072973525643

def nonUniqueShapeFactor(Z,A,Q_kev,Ji,Jf,chemin_test):

    Q = Q_kev/ (Me*1000)
    W0 = 1 + Q
    R = 3.11e-3 * A**(1/3)
    
    n = int(NUM_BINS * Q / (E_MAX/Me))
    W = np.linspace(1, W0, n)
    p = np.sqrt(W**2.0-1.)  # confirmed from thecobs/src/coulomb_functions.py
    q = W0 - W              # momentum of neutrino p_nu = W_nu = q from Behrens Buhring

    def FA(N, K, L, s):
        return FANKLs(N, K, L, s, 0,0,0,0, chemin_test)
        
    def FV(N, K, L, s):
        return FVNKLs(N, K, L, s, 0,0,0,0, chemin_test)

    def FA_1111(N, K, L, s):
        return FANKLs(N, K, L, s, 1,1,1,1, chemin_test)
        
    def FV_1111(N, K, L, s):
        return FVNKLs(N, K, L, s, 1,1,1,1, chemin_test)



    def M0_11():
        return FA(0,0,0,0) - 1/3*ALPHA*Z*FA_1111(0,0,1,1) - 1/3*W0*R*FA(0,0,1,1)

    def m0_11():
        return -1/3*(1)*R*FA(0,0,1,1) # m_e = 1 in natural units

    def M1_11():
        return -FV(0,1,0,1) - 1/3*ALPHA*Z*(1/3)**(0.5)*FV_1111(0,1,1,0) \
        - 1/3*W0*R*(1/3)**(0.5)*FV(0,1,1,0) - 1/3*ALPHA*Z*(2/3)**(0.5)*FA_1111(0,1,1,1) \
        - 1/3*(W-q)*R*(2/3)**(0.5)*FA(0,1,1,1) 

    def m1_11():
        return -1/3*(1)*R*((1/3)**(0.5)*FV(0,1,1,0) + (2/3)**(0.5)* FA(0,1,1,1)) # m_e = 1 in natural units

    def M1_12():
        return 1/3*q*R*((2/3)**(0.5)*FV(0,1,1,0) + (1/3)**(0.5)*FA(0,1,1,1))

    def M1_21():
        return 1/3*p*R*((2/3)**(0.5)*FV(0,1,1,0) - (1/3)**(0.5)*FA(0,1,1,1))

    def M2_12():
        return -1/3*q*R*FA(0,2,1,1)
        
    def M2_21():
        return -1/3*p*R*FA(0,2,1,1)



    def m(L,k_e,k_nu):
        # print("m(",L,k_e,k_nu,")")
        if L == 0 and k_e == 1 and k_nu == 1:
            return m0_11()
        elif L == 1 and k_e == 1 and k_nu == 1:
            return m1_11()
        # else:
        #     print("neglecting: m(",L,k_e,k_nu,")",sep=",")
        return 0

    def M(L,k_e,k_nu):
        # print("M(",L,k_e,k_nu,")")
        if L == 0 and k_e == 1 and k_nu == 1:
            return M0_11()
        
        elif L == 1 and k_e == 1 and k_nu == 1:
            return M1_11()
        elif L == 1 and k_e == 1 and k_nu == 2:
            return M1_12()
        elif L == 1 and k_e == 2 and k_nu == 1:
            return M1_21()
        
        elif L == 2 and k_e == 1 and k_nu == 2:
            return M2_12()
        elif L == 2 and k_e == 2 and k_nu == 1:
            return M2_21()
        # else:
        #     print("neglecting: M(",L,k_e,k_nu,")",sep=",")
        return 0

    def C():
        if Jf - Ji == 0:
            L = 1
        else:
            L = abs(Jf-Ji)
        
        term1 = 0.
        for k_e in range(1,int(L)+1):
            k_nu = L+1 - k_e

            lambda_ke = Utilities.lambda_ke(k_e, W, chemin_test)
            mu_ke = Utilities.mu_ke(k_e, W, chemin_test)

            gamma_ke = (k_e**2 - (ALPHA*Z)**2)**(0.5)

            M_L_ke_knu = M(L,k_e,k_nu)
            m_L_ke_knu = m(L,k_e,k_nu)

            mM_term = 2*mu_ke*gamma_ke/k_e/W * M_L_ke_knu * m_L_ke_knu

            term1 += lambda_ke * (M_L_ke_knu**2 + m_L_ke_knu**2 - mM_term)
        
        term2 = 0.
        
        for k_e in range(1,int(L)+2):
            k_nu = L+2 - k_e

            lambda_ke = Utilities.lambda_ke(k_e, W, chemin_test)
            mu_ke = Utilities.mu_ke(k_e, W, chemin_test)

            gamma_ke = (k_e**2 - (ALPHA*Z)**2)**(0.5)

            M_L_ke_knu = M(L,k_e,k_nu)
            M_L1_ke_knu = M(L+1,k_e,k_nu)
            m_L_ke_knu = m(L,k_e,k_nu)
            m_L1_ke_knu = m(L+1,k_e,k_nu)

            mM_term = 2*mu_ke*gamma_ke/k_e/W * M_L_ke_knu * m_L_ke_knu
            mM_term2 = 2*mu_ke*gamma_ke/k_e/W * M_L1_ke_knu * m_L1_ke_knu

            term2 += lambda_ke*(M_L_ke_knu**2 - mM_term + M_L1_ke_knu**2 - mM_term2)
        
        if Jf - Ji == 0:
            k_e = 1
            mu_ke = Utilities.mu_ke(1, W, chemin_test)
            gamma_ke = (k_e**2 - (ALPHA*Z)**2)**(0.5)
            M011 = M(0,1,1)
            m011 = m(0,1,1)
            return term1 + term2 + M011**2 + m011**2 - 2*mu_ke*gamma_ke/k_e/W * M011 * m011
        else:
            return term1 + term2
        
    return W, C()

# W, shape_factor = nonUniqueShapeFactor(83, 210, 1162.2, 1.0, 0.0, "./OBTDs/test-Bi210.txt")

# import matplotlib.pyplot as plt

# plt.figure('')
# plt.plot(W, shape_factor)
# plt.xlabel("W (Natural Units)")
# plt.ylabel("C(E_e)")
# plt.title(f"shape factor")
# plt.legend()
# plt.grid()
# plt.show()