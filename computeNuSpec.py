import numpy as np
import thecobs.SpectralFunctions as sf
import radioactivedecay as rd

from Resolution import *
from shapeFactor import nonUniqueShapeFactor
from Helper import getJiJf

np.seterr(divide='ignore', invalid='ignore')

def betaSpectrum(Z,A,Q): # everything happens here in natural units
    W0 = 1 + Q # already converted from keV to electron rest mass units

    n = int(NUM_BINS * Q / (E_MAX/Me))
    W0 = n * (E_MAX/Me) / NUM_BINS + 1

    R = 3.11e-3 * A**(1/3)         # in electron compton units
    W = np.linspace(1, W0, n)
    l = 2/137*Z**0.33

    ps = sf.phase_space(W, W0)
    f = sf.fermi_function(W, Z, R)
    g = sf.sirlin_g(W, W0)
    rf = sf.recoil_fermi(W, W0, A)
    rc1 = (1+sf.ALPHA/2/np.pi*g)
    rc = sf.radiative_correction(W, Z, W0, R)
    s = sf.atomic_screening(W, Z, R, l)

    E = (W-1)
    intensity = ps*f*rc*rf*s
    return intensity, W, W0, Z, R


def computeNuSpec(Transitions, shapeFactorNU = {}, random = False):

    spectra = {}
    for pr_Nuc in Transitions.keys():
        beta = dict()
        beta = Transitions[pr_Nuc]

        pr_spec = np.zeros(shape=(1))

        for nuc, peaks in beta.items():

            dn_spec = np.zeros(shape=(1))

            for Q_kev, info in peaks.items():

                br_ratio = info['Contribution']

                # args in natural units
                A = rd.Nuclide(nuc).A
                intensity, W, W0, Z, R = betaSpectrum(
                                            rd.Nuclide(nuc).Z,
                                            rd.Nuclide(nuc).A,
                                            Q_kev/ (Me*1000)
                                        )

                ## apply shape factors

                # intensity = intensity * sf.shape_factor_fermi(W, Z, W0, R)

                if info['Type'] == 'A':
                  C_w = sf.shape_factor_gamow_teller(W, Z, W0, R, A, 4.5*A, 1,0, 0)
                  intensity = intensity * C_w

                  # ## Plot shape factors
                  # plt.plot(E[0:len(sf_npy)], sf_npy)
                  # plt.title(nuc + f" Allowed ({info['j_ini']} to {info['j_fin']}), Endpoint: " + str(round(Q_kev/1000,2)) + " MeV")
                  # plt.xlabel('Electron Energy (MeV)')
                  # plt.ylabel('Shape Factor C(E)')
                  # plt.show()

                elif info['Type'][1:] == 'NU':

                    Ji,Jf = getJiJf(info['Ji'], info['Jf'])

                    if Ji is None or Jf is None:
                        pass

                    elif (nuc,Z,A,Q_kev,Ji,Jf) in shapeFactorNU.keys():
                        print("Computing NU forbidden shape factor:", (nuc,Z,A,Q_kev,Ji,Jf))
                        chemin_test = shapeFactorNU[(nuc,Z,A,Q_kev,Ji,Jf)]

                        _, C_w = nonUniqueShapeFactor(Z, A, Q_kev, Ji, Jf, chemin_test)
                        C_w[0] = C_w[1]

                        extra_ones = len(intensity) - len(C_w)
                        C_w = np.concatenate((C_w, np.array([C_w[-1]]*extra_ones)))

                        C_w = C_w[0:len(intensity)]
                        intensity = intensity * C_w


                    elif info['Type'] == '1NU' and random:
                      
                        cmat = [[1.00, -0.62,  -0.98],
                                [-0.62, 1.00,  0.55],
                                [-0.98,  0.55,  1.00]]

                        # print(info['j_ini'], info['j_fin'])

                        L = abs(Ji - Jf)
                        do = True
                        if L == 1: # orange
                            mean = np.array([-0.096, 0.14, 0.006])
                            std = np.array([0.13, 0.18, 0.014])

                        elif L == 0: # blue
                            mean = np.array([-0.005, 0.0, 0.001])
                            std = np.array([0.050, 0.054, 0.0045])
                        else:
                            do = False

                        if do :
                            sigma = np.diag(std)
                            cov_mat = np.matmul(np.matmul(sigma, cmat), sigma)

                            coeff = np.random.multivariate_normal(mean, cov_mat)
                            a = coeff[0]
                            b = coeff[1]
                            c = coeff[2]
                            intensity = intensity * (1 + a*W + b/W + c*np.square(W))


                elif info['Type'][1:] == 'U':
                    if not (info['Ji'] == '' or info['Jf'] == ''):

                        Ji,Jf = getJiJf(info['Ji'], info['Jf'])
                        L = int(abs(Ji - Jf))
                        C_w = sf.shape_factor_unique_forbidden(W, L, W0, Z, R)
                        intensity = intensity * C_w

                        # plt.plot(E[0:len(sf_npy)], sf_npy)
                        # plt.title(nuc + f" {info['type']} ({info['j_ini']} to {info['j_fin']}), Endpoint: " + str(round(Q_kev/1000,2)) + " MeV")
                        # plt.xlabel('Electron Energy (MeV)')
                        # plt.ylabel('Shape Factor C(E)')
                        # plt.show()


                intensity[np.isnan(intensity)] = 0.
                intensity = np.flip(intensity, axis=0) # to get neutrino spectra

                extra_zeros = NUM_BINS - len(intensity)
                intensity = np.concatenate((intensity, np.zeros(shape=extra_zeros)))

                # normalized spectrum
                norm = np.sum(intensity) * dE
                norm_spec = intensity / norm

                dn_spec = dn_spec + (br_ratio * norm_spec)

            if len(dn_spec) > 1:
              pr_spec = pr_spec + dn_spec

        spectra[pr_Nuc] = pr_spec

    return E, spectra
