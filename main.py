from computeNuSpec import computeNuSpec
from queryDecayChain import queryDecayChain

transitions = dict()
nuclei = ['U-238', 'Th-232', 'U-235',  'K-40', 'Rb-87']

# for pr_Nuc in nuclei:
#     transitions[pr_Nuc] = dict()
#     transitions[pr_Nuc] = queryDecayChain(pr_Nuc, 1., transitions[pr_Nuc])

from transitions import transitions

# Don't know why not to include the transition: Bi-212  807.0   1NU     (18-)   (18+)
transitions["Th-232"]["Bi-212"][807.0]["Contribution"] = 0.

# ## Print all Transitions  
# for prNuc in transitions.keys():
#     betaNuclei = transitions[prNuc]
#     transitionData = []
#     areaCurve = 0.
#     for dNuc, decayData in betaNuclei.items():

#         for Qval, data in decayData.items():

#             areaCurve += data['Contribution']
            
#             percentage = round(data['Contribution'] * 100, 3)
#             transitionData.append((percentage, dNuc, Qval, data['Type'], data['Ji'], data['Jf']))

#     transitionData.sort()
#     transitionData.reverse()

#     print(prNuc, round(areaCurve * 100,2), '%')
#     for peak in transitionData:

#         print('\t'.join([str(i) for i in peak]))
#     print('##################\n\n')

shapeFactorNU = {}

shapeFactorNU[("Bi-210", 83, 210, 1162.2, 1.0, 0.0)] = "./OBTDs/test-Bi210.txt"
shapeFactorNU[("Bi-212", 83, 212, 2251.5, 1.0, 0.0)] = "./OBTDs/test-Bi212.txt"
shapeFactorNU[("Bi-214", 83, 214, 3269.0, 1.0, 0.0)] = "./OBTDs/test-Bi214.txt"
shapeFactorNU[("Bi-215", 83, 215, 2189.0, 9/2, 9/2)] = "./OBTDs/test-Bi215.txt"

# ## Compute Spectra
E, spectra = computeNuSpec(transitions, shapeFactorNU=shapeFactorNU)


# ## Plot Spectra
import matplotlib.pyplot as plt

colors = ['blue', 'red', 'green', 'magenta', 'yellow']
for i in range(len(nuclei)):
    N_nu = spectra[nuclei[i]]
    plt.plot(E, N_nu, color = colors[i])

# ## Threshold 1.8 MeV
plt.plot([1.8,1.8], [0.0005, 300], linestyle='dotted', color='black')

plt.legend(nuclei)
plt.yscale('log')
plt.xlabel('Neutrino Energy (MeV)')
plt.ylabel('Intensity (decay^-1 MeV^-1)')

plt.xlim(right = 3.5)
plt.ylim(top = 300)  # adjust the top leaving bottom unchanged
plt.ylim(bottom = 0.0005)  # adjust the bottom leaving top unchanged
plt.show()
