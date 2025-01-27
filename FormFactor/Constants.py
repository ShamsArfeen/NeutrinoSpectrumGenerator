#################################################
# Constants.py file
#################################################
import numpy as np

#################################################
# Physic Constants
#################################################
me = 510.9989461 #keV, electron mass
mp = 938272.046 #keV, proton mass
mn = 939565.4133 #keV, neutron mass
UMASSC2 = 931494.028 #kev
ALPHA = 1/137 # fine structure constante
e = np.sqrt( (4 * np.pi * ALPHA))
charge_e = 1.6021766208 * 10**(-19) # C
amm_p = 2.79284734463 # anomalous magnetic moment of proton
amm_n = -1.91304276 # # anomalous magnetic moment of proton
#################################################

#################################################
# Physic Parameters
#################################################
GF = 1.166378 * 10**(-5)
Vud = 0.97435
Gv = GF*Vud

q = 0.38 # quenching for 210-215Bi, in according with Shell-model description for the properties of the forbidden β− decay in the region “northeast” of 208Pb, Shweta Sharma, Praveen C. Srivastava, Anil Kumar and Toshio Suzuki
gv = 1 
LAMBDA = ga = -1.27  * q # Axial coupling * quenching


fM = 0 #(amm_p-1-amm_n)/(2*UMASSC2/me)
fS = 0
fT = 0
fP = 0
#################################################

#################################################
# Nuclear occupation according with BSG code
#################################################
sm_occupation = [
    1, 0, 1, 2,   1, 1, 1, 6,   1, 1, -1, 8,   1, 2, 1, 14,   2, 0, 1, 16,
    1, 2, -1, 20,  1, 3, 1, 28,   2, 1, 1, 32,  1, 3, -1, 38,   2, 1, -1, 40,
    1, 4, 1, 50,   1, 4, -1, 58,  2, 2, 1, 64,   2, 2, -1, 68,   3, 0, 1, 70,
    1, 5, 1, 82,   1, 5, -1, 92,  2, 3, 1, 100,  2, 3, -1, 106,  3, 1, 1, 110,
    3, 1, -1, 112,  1, 6, 1, 126,  2, 4, 1, 136,  3, 2, 1, 142,  1, 6, -1, 154
]
#################################################

#################################################
# For 210Bi case
#################################################
E0 = 3.27
Z = 83

