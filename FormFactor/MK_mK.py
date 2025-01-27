#################################################
# MK_mK.py file 
#################################################
import numpy as np

from scipy.special import factorial, factorial2


#################################################
# Kn factor in according with Nuclear Beta Decay, 1971, H.Behrens and W.Bühring
#################################################
def Kn(n, ke, knu):
    kn = np.sqrt(1/2) * np.sqrt(factorial2(2*n)/factorial2(2*n+1)) * np.sqrt(1/(factorial(2*ke-1)*factorial(2*knu-1)))
    return kn
