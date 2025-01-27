#################################################
# FF.py file (Form Factor)
#################################################
import numpy as np

import Utilities
import ME
import IF
import NWF
import Constants

E0 = Constants.E0
#################################################
# Form Factor
#################################################
def FVNKK0(N, K, L, s, k, m, n, rho, chemin) :
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 30*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    mv, mm, ms = 0, 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '0') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)

                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1

                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_m12 = []
                psi_m34 = []
                psi_s12 = []
                for i in range(len(r)):
                    psi_m12.append((2*K+1+2*N)*II[i]+r[i]*I_prime[i])
                    psi_m34.append(2*N*II[i]+r[i]*I_prime[i])
                    psi_s12.append(E0*R/385-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctv1 = gf*pow(r/R, K+2*N)*II*gi*pow(r, 2)
                fctv2 = ff*pow(r/R, K+2*N)*II*fi*pow(r, 2)

                fctm1 = gf*pow(r/R, K-1+2*N)*psi_m12*fi*pow(r, 2)
                fctm2 = ff*pow(r/R, K-1+2*N)*psi_m12*gi*pow(r, 2)
                fctm3 = gf*pow(r/R, K-1+2*N)*psi_m34*fi*pow(r, 2)
                fctm4 = ff*pow(r/R, K-1+2*N)*psi_m34*gi*pow(r, 2)

                fcts1 = gf*pow(r/R, K+2*N)*II*psi_s12*gi*pow(r, 2)
                fcts2 = ff*pow(r/R, K+2*N)*II*psi_s12*fi*pow(r, 2)

                intv1 = 0
                intv2 = 0
                intm1, intm2, intm3, intm4 = 0, 0, 0, 0
                ints1, ints2 = 0, 0
                for i in range(len(r)-1):
                    intv1 += (fctv1[i]+fctv1[i+1])/2 * dr
                    intv2 += (fctv2[i]+fctv2[i+1])/2 * dr

                    intm1 += (fctm1[i]+fctm1[i+1])/2 * dr
                    intm2 += (fctm2[i]+fctm2[i+1])/2 * dr
                    intm3 += (fctm3[i]+fctm3[i+1])/2 * dr
                    intm4 += (fctm4[i]+fctm4[i+1])/2 * dr

                    ints1 += (fcts1[i]+fcts1[i+1])/2 * dr
                    ints2 += (fcts2[i]+fcts2[i+1])/2 * dr
                
                mv += ME.MVNKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mm += Constants.fM/R * (np.sqrt(K/(2*K+1))*np.sqrt(2)/np.sqrt(2*Ji+1) * (-np.sign(ki)*Utilities.GKLs(K, L-1, 1, kf, -ki)*intm1 + np.sign(kf)*Utilities.GKLs(K, L-1, 1, -kf, ki)*intm2) + np.sqrt((K+1)/(2*K+1))*np.sqrt(2)/np.sqrt(2*Ji+1) * (-np.sign(ki)*Utilities.GKLs(K, L+1, 1, kf, -ki)*intm3 + np.sign(kf)*Utilities.GKLs(K, L+1, 1, -kf, ki)*intm4)) * obd
                ms += Constants.fS/R * np.sqrt(2)/np.sqrt(2*Ji+1)*(-Utilities.GKLs(K, L, s, kf, ki)*ints1 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, L, s, -kf, -ki)*ints2) * obd

        f.close()
    if beta_type == 'beta+' :
        fvnkk0 = float(mv) + float(mm) + float(ms)
    else :
        fvnkk0 =float(mv) + float(mm) - float(ms)
    return fvnkk0 

def FVNKK1(N, K, L, s, k, m, n, rho, chemin):
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 30*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    mv, mm = 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '0') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)
                
                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1
        
                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_m12 = []
                psi_m34 = []
                psi_m56 = []
                for i in range(len(r)):
                    psi_m12.append((2*K+1+2*N)*II[i]+r[i]*I_prime[i])
                    psi_m34.append(2*N*II[i]+r[i]*I_prime[i])
                    psi_m56.append(E0*R-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctv1 = gf*pow(r/R, K+2*N)*II*fi*pow(r, 2)
                fctv2 = ff*pow(r/R, K+2*N)*II*gi*pow(r, 2)

                fctm1 = gf*pow(r/R, K-1+2*N)*psi_m12*gi*pow(r, 2)
                fctm2 = ff*pow(r/R, K-1+2*N)*psi_m12*fi*pow(r, 2)
                fctm3 = gf*pow(r/R, K-1+2*N)*psi_m34*gi*pow(r, 2)
                fctm4 = ff*pow(r/R, K-1+2*N)*psi_m34*fi*pow(r, 2)

                fctm5 = gf*pow(r/R, L+2*N)*II*psi_m56*fi*pow(r, 2)
                fctm6 = ff*pow(r/R, L+2*N)*II*psi_m56*gi*pow(r, 2)

                intv1, intv2 = 0, 0
                intm1, intm2, intm3, intm4 = 0, 0, 0, 0
                intm5, intm6 = 0, 0
                for i in range(len(r)-1):
                    intv1 += (fctv1[i]+fctv1[i+1])/2 * dr
                    intv2 += (fctv2[i]+fctv2[i+1])/2 * dr

                    intm1 += (fctm1[i]+fctm1[i+1])/2 * dr
                    intm2 += (fctm2[i]+fctm2[i+1])/2 * dr
                    intm3 += (fctm3[i]+fctm3[i+1])/2 * dr
                    intm4 += (fctm4[i]+fctm4[i+1])/2 * dr

                    intm5 += (fctm5[i]+fctm5[i+1])/2 * dr
                    intm6 += (fctm6[i]+fctm6[i+1])/2 * dr
                
                mv += ME.MVNKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mm += Constants.fM/R*np.sqrt(2)/np.sqrt(2*Ji+1)* ( (np.sqrt(K/(2*K+1))*(-Utilities.GKLs(K, K-1, 1, kf, ki)*intm1 + np.sign(ki)*np.sign(kf) * Utilities.GKLs(K, K-1, 1, -kf, -ki)*intm2) - np.sqrt((K+1)/(2*K+1))*(-Utilities.GKLs(K, K+1, 1, kf, ki)*intm3 + np.sign(ki)*np.sign(kf) * Utilities.GKLs(K, K+1, 1, -kf, -ki)*intm4)) + (-np.sign(ki)*Utilities.GKLs(K, K, 1, kf, -ki)*intm5 + np.sign(kf)*Utilities.GKLs(K, K, 1, -kf, ki)*intm6) ) * obd
       
        f.close()    
    if beta_type == 'beta+' :
        fvnkk1 = float(mv) + float(mm)
    else :
        fvnkk1 =float(mv) + float(mm)
    return fvnkk1

def FVNKLm1(N, K, L, s, k, m, n, rho, chemin):
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 30*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    mv, mm, ms = 0, 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '0') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)
                
                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1
        
                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_m12 = []
                psi_m34 = []
                for i in range(len(r)):
                    psi_m12.append(2*N*II[i]+r[i]*I_prime[i])
                    psi_m34.append(E0*R-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctm1 = gf*pow(r/R, K-3+2*N)*psi_m12*gi*pow(r, 2)
                fctm2 = ff*pow(r/R, K-3+2*N)*psi_m12*fi*pow(r, 2)
                fctm3 = gf*pow(r/R, K-1+2*N)*II*psi_m34*gi*pow(r, 2)
                fctm4 = ff*pow(r/R, K-1+2*N)*II*psi_m34*fi*pow(r, 2)

                fcts1 = gf*pow(r/R, K-2+2*N)*psi_m12*gi*pow(r, 2)
                fcts2 = ff*pow(r/R, K-2+2*N)*psi_m12*fi*pow(r, 2)

                intm1, intm2, intm3, intm4 = 0, 0, 0, 0
                ints1, ints2 = 0, 0
                for i in range(len(r)-1):
                    intm1 += (fctm1[i]+fctm1[i+1])/2 * dr
                    intm2 += (fctm2[i]+fctm2[i+1])/2 * dr
                    intm3 += (fctm3[i]+fctm3[i+1])/2 * dr
                    intm4 += (fctm4[i]+fctm4[i+1])/2 * dr

                    ints1 += (fcts1[i]+fcts1[i+1])/2 * dr
                    ints2 += (fcts2[i]+fcts2[i+1])/2 * dr
                    
                mv += ME.MVNKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mm += Constants.fM/R*np.sqrt(2)/np.sqrt(2*Ji+1)* (-np.sqrt((K+1)/(2*K+1))*(-Utilities.GKLs(K, K, 1, kf, ki)*intm1 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, K, 1, -kf, -ki)*intm2) + (-Utilities.GKLs(K, K-1, 1, kf, ki)*intm3 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, K-1, 1, -kf, -ki)*intm4)) * obd
                ms += -Constants.fS/R*np.sqrt(K/(2*K+1))*np.sqrt(2)/np.sqrt(2*Ji+1)*(-Utilities.GKLs(K, K, 0, kf, ki)*ints1 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, K, 0, -kf, -ki)*ints2)

        f.close()
    if beta_type == 'beta+' :
        fvnklm1 = -(float(mv) + float(mm) + float(ms))
    else :
        fvnklm1 = -(float(mv) + float(mm) - float(ms))
    return fvnklm1

def FVNKLp1(N, K, L, s, k, m, n, rho, chemin):
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 30*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    mv, mm, ms = 0, 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '0') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)
                
                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1
        
                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_m12 = []
                psi_m34 = []
                for i in range(len(r)):
                    psi_m12.append((2*K+3+2*N)*II[i]+r[i]*I_prime[i])
                    psi_m34.append(E0*R-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctm1 = gf*pow(r/R, K+2*N)*psi_m12*gi*pow(r, 2)
                fctm2 = ff*pow(r/R, K+2*N)*psi_m12*fi*pow(r, 2)
                fctm3 = gf*pow(r/R, K+1+2*N)*II*psi_m34*fi*pow(r, 2)
                fctm4 = ff*pow(r/R, K+1+2*N)*II*psi_m34*gi*pow(r, 2)

                fcts1 = gf*pow(r/R, K+2*N)*psi_m12*gi*pow(r, 2)
                fcts2 = ff*pow(r/R, K+2*N)*psi_m12*fi*pow(r, 2)

                intm1, intm2, intm3, intm4 = 0, 0, 0, 0
                ints1, ints2 = 0, 0
                dr = r[1] - r[0]
                for i in range(len(r)-1):
                    intm1 += (fctm1[i]+fctm1[i+1])/2 * dr
                    intm2 += (fctm2[i]+fctm2[i+1])/2 * dr
                    intm3 += (fctm3[i]+fctm3[i+1])/2 * dr
                    intm4 += (fctm4[i]+fctm4[i+1])/2 * dr

                    ints1 += (fcts1[i]+fcts1[i+1])/2 * dr
                    ints2 += (fcts2[i]+fcts2[i+1])/2 * dr
                
                mv += ME.MVNKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mm += Constants.fM/R*np.sqrt(2)/np.sqrt(2*Ji+1)* (np.sqrt(K/(2*K+1))*(-Utilities.GKLs(K, K, 1, kf, ki)*intm1 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, K, 1, -kf, -ki)*intm2 ) + (-np.sign(ki)*Utilities.GKLs(K, K+1, 1, kf, -ki)*intm3 + np.sign(kf)*Utilities.GKLs(K, K+1, 1, -kf, ki)*intm4)) * obd
                ms += -Constants.fS/R*np.sqrt((K+1)/(2*K+1))*np.sqrt(2)/np.sqrt(2*Ji+1)* (-Utilities.GKLs(K, K, 0, kf, ki)*ints1 + np.sign(kf)*np.sign(ki)*Utilities.GKLs(K, K, 0, -kf, -ki)*ints2) * obd

        f.close()
    if beta_type == 'beta+' :
        fvnklp1 = -(float(mv) + float(mm) + float(ms))
    else :
        fvnklp1 = -(float(mv) + float(mm) - float(ms))
    return fvnklp1

def FANKK0(N, K, L, s, k, m, n, rho, chemin):
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 30*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    ma, mt, mp = 0, 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '1') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)
                
                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1
        
                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_t12 = []
                psi_t34 = []
                psi_p12 = []
                for i in range(len(r)):
                    psi_t12.append((2*K+1+2*N)*II[i]+r[i]*I_prime[i])
                    psi_t34.append(2*N*II[i]+r[i]*I_prime[i])
                    psi_p12.append(E0*R-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctt1 = gf*pow(r/R, K-1+2*N)*psi_t12*gi*pow(r, 2)
                fctt2 = ff*pow(r/R, K-1+2*N)*psi_t12*fi*pow(r, 2)
                fctt3 = gf*pow(r/R, K-1+2*N)*psi_t34*gi*pow(r, 2)
                fctt4 = ff*pow(r/R, K-1+2*N)*psi_t34*fi*pow(r, 2)

                fctp1 = gf*pow(r/R, K+2*N)*II*psi_p12*fi*pow(r, 2)
                fctp2 = ff*pow(r/R, K+2*N)*II*psi_p12*gi*pow(r, 2)

                intt1, intt2, intt3, intt4 = 0, 0, 0, 0
                intp1, intp2 = 0, 0
                dr = r[1] - r[0]
                for i in range(len(r)-1):
                    intt1 += (fctt1[i]+fctt1[i+1])/2 * dr
                    intt2 += (fctt2[i]+fctt2[i+1])/2 * dr
                    intt3 += (fctt3[i]+fctt3[i+1])/2 * dr
                    intt4 += (fctt4[i]+fctt4[i+1])/2 * dr

                    intp1 += (fctp1[i]+fctp1[i+1])/2 * dr
                    intp2 += (fctp2[i]+fctp2[i+1])/2 * dr
                
                ma += ME.MANKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mt += - Constants.fT/R * (np.sqrt(K/(2*K+1))*np.sqrt(2)/np.sqrt(2*Ji+1)*(-Utilities.GKLs(K, K-1, 1, kf, ki)*intt1 + np.sign(ki)*np.sign(kf) * Utilities.GKLs(K, K-1, 1, -kf, -ki)*intt2) + np.sqrt((K+1)/(2*K+1))*np.sqrt(2)/np.sqrt(2*Ji+1)*(-Utilities.GKLs(K, K+1, 1, kf, ki)*intt3 + np.sign(ki)*np.sign(kf) * Utilities.GKLs(K, K+1, 1, -kf, -ki)*intt4) ) * obd
                mp += - Constants.fP/R * np.sqrt(2)/np.sqrt(2*Ji+1) * (-np.sign(ki)*Utilities.GKLs(K, L, 0, kf, -ki)*intp1 + np.sign(kf)*Utilities.GKLs(K, L, 0, -kf, ki)*intp2) * obd

        f.close()
    if beta_type == 'beta+' :
        fankk0 = Constants.LAMBDA * float(ma) + float(mt) + float(mp)
    else :
        fankk0 = Constants.LAMBDA * float(ma) - float(mt) + float(mp)
    return fankk0

def FANKK1(N, K, L, s, k, m, n, rho, chemin):
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 15*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    ma, mt, mp = 0, 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '0') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)

                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1

                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_t12 = []
                psi_t34 = []
                psi_t56 = []
                for i in range(len(r)):
                    psi_t12.append((2*K+1+2*N)*II[i]+r[i]*I_prime[i])
                    psi_t34.append(2*N*II[i]+r[i]*I_prime[i])
                    psi_t56.append(E0*R/385-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctt1 = gf*pow(r/R, K-1+2*N)*psi_t12*fi*pow(r, 2)
                fctt2 = ff*pow(r/R, K-1+2*N)*psi_t12*gi*pow(r, 2)
                fctt3 = gf*pow(r/R, K-1+2*N)*psi_t34*fi*pow(r, 2)
                fctt4 = ff*pow(r/R, K-1+2*N)*psi_t34*gi*pow(r, 2)
                fctt5 = gf*pow(r/R, K+2*N)*II*psi_t56*gi*pow(r, 2)
                fctt6 = ff*pow(r/R, K+2*N)*II*psi_t56*fi*pow(r, 2)

                intt1, intt2, intt3, intt4, intt5, intt6 = 0, 0, 0, 0, 0, 0
                dr = r[1] - r[0]
                for i in range(len(r)-1):
                    intt1 += (fctt1[i]+fctt1[i+1])/2 * dr
                    intt2 += (fctt2[i]+fctt2[i+1])/2 * dr
                    intt3 += (fctt3[i]+fctt3[i+1])/2 * dr
                    intt4 += (fctt4[i]+fctt4[i+1])/2 * dr
                    intt5 += (fctt5[i]+fctt5[i+1])/2 * dr
                    intt6 += (fctt6[i]+fctt6[i+1])/2 * dr

                ma += ME.MANKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mt += -Constants.fT/R*np.sqrt(2)/np.sqrt(2*Ji+1)* (np.sqrt((K+1)/(2*K+1))*(-np.sign(ki)*Utilities.GKLs(K, K-1, 1, kf, -ki)*intt1 + np.sign(kf)*Utilities.GKLs(K, K-1, 1, -kf, ki)*intt2) - np.sqrt(K/(2*K+1))*(-np.sign(ki)*Utilities.GKLs(K, K+1, 1, kf, -ki)*intt3 + np.sign(kf)*Utilities.GKLs(K, K+1, 1, -kf, ki)*intt4) + (-Utilities.GKLs(K, K, 1, kf, ki)*intt5 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, K, 1, -kf, -ki)*intt6)) * obd

        f.close()
    if beta_type == 'beta+' :
        fankk1 = Constants.LAMBDA * float(ma) + float(mt)
    else :
        fankk1 = Constants.LAMBDA * float(ma) - float(mt)
    return fankk1
    
def FANKLm1(N, K, L, s, k, m, n, rho, chemin):
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 30*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    ma, mt, mp = 0, 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '1') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)
                
                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1
        
                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_t12 = []
                psi_t34 = []
                for i in range(len(r)):
                    psi_t12.append(2*N*II[i]+r[i]*I_prime[i])
                    psi_t34.append(E0*R-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctt1 = gf*pow(r/R, K-2+2*N)*psi_t12*fi*pow(r, 2)
                fctt2 = ff*pow(r/R, K-2+2*N)*psi_t12*gi*pow(r, 2)
                fctt3 = gf*pow(r/R, K-1+2*N)*II*psi_t34*gi*pow(r, 2)
                fctt4 = ff*pow(r/R, K-1+2*N)*II*psi_t34*fi*pow(r, 2)
                fctp1 = gf*pow(r/R, K-2+2*N)*psi_t12*fi*pow(r, 2)
                fctp2 = ff*pow(r/R, K-2+2*N)*psi_t12*gi*pow(r, 2)

                intt1, intt2, intt3, intt4 = 0, 0, 0, 0
                intp1, intp2 = 0, 0
                dr = r[1] - r[0]
                for i in range(len(r)-1):
                    intt1 += (fctt1[i]+fctt1[i+1])/2 * dr
                    intt2 += (fctt2[i]+fctt2[i+1])/2 * dr
                    intt3 += (fctt3[i]+fctt3[i+1])/2 * dr
                    intt4 += (fctt4[i]+fctt4[i+1])/2 * dr
                    intp1 += (fctp1[i]+fctp1[i+1])/2 * dr
                    intp2 += (fctp2[i]+fctp2[i+1])/2 * dr
                
                ma += ME.MANKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mt += -Constants.fT/R*np.sqrt(2)/np.sqrt(2*Ji+1)* (-np.sqrt((K+1)/(2*K+1))*(-np.sign(ki)*Utilities.GKLs(K, K, 1, kf, -ki)*intt1 + np.sign(kf)*Utilities.GKLs(K, K, 1, -kf, ki)*intt2) + (-Utilities.GKLs(K, K-1, 1, kf, ki)*intt3 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, K-1, 1, -kf, -ki)*intt4) ) * obd
                mp += Constants.fP/R*np.sqrt(2)/np.sqrt(2*Ji+1)*np.sqrt(K/(2*K+1)) * (-np.sign(ki)*Utilities.GKLs(K, K, 0, kf, -ki)*intp1 + np.sign(kf)*Utilities.GKLs(K, K, 0, -kf, ki)*intp2) * obd

        f.close()
    if beta_type == 'beta+' :
        fanklm1 = -(Constants.LAMBDA * float(ma) + float(mt) + float(mp))
    else :
        fanklm1 = -(Constants.LAMBDA * float(ma) - float(mt) + float(mp))
    return fanklm1
    
def FANKLp1(N, K, L, s, k, m, n, rho, chemin):
    beta_type, A, Z_mother, Ji, nucleus_name = Utilities.find_nucleus(chemin)
    R = 1.2 * A**(1/3) # (1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A)
    r = np.linspace(1e-6, 30*R, 5000)
    dr = r[1] - r[0]
    if beta_type == 'beta+' :
        Z_daughter = Z_mother - 1
    else :
        Z_daughter = Z_mother + 1
        
    ma, mt, mp = 0, 0, 0
    with open(chemin, newline='\n') as fichier_txt: 
        f = open(chemin, 'r')
        for ligne in fichier_txt :
            data = ligne.split()

            if (len(data) == 6) : # and (data[0] == '1') :
                ki, ni = int(data[1]), int(data[3])
                kf, nf = int(data[2]), int(data[4])
                obd = float(data[5])
                
                nui, nuf = NWF.calc_nu(R, Z_mother), NWF.calc_nu(R, Z_daughter)
                
                lnf = kf
                lni = ki

                if lnf < 0 :
                    lnf = abs(lnf) - 1
                if lni < 0 : 
                    lni = abs(lni) - 1
        
                gf, ff = NWF.NWF(nf, lnf, nuf, r)
                gi, fi = NWF.NWF(ni, lni, nui, r)
                
                II = []
                I_prime= []
                for i in range(len(r)):
                    II.append(IF.I(k, r[i], m, n, rho, R))
                for i in range(len(r)-1):
                    I_prime.append((II[i+1]-II[i])/dr)
                I_prime.append(I_prime[-1])

                psi_t12 = []
                psi_t34 = []
                for i in range(len(r)):
                    psi_t12.append((2*K+3+2*N)*II[i]+r[i]*I_prime[i])
                    psi_t34.append(E0*R-Constants.ALPHA*Z_mother*Utilities.Ur(r[i], Z_mother, R))

                fctt1 = gf*pow(r/R, K+2*N)*psi_t12*fi*pow(r, 2)
                fctt2 = ff*pow(r/R, K+2*N)*psi_t12*gi*pow(r, 2)
                fctt3 = gf*pow(r/R, K+1+2*N)*II*psi_t34*gi*pow(r, 2)
                fctt4 = ff*pow(r/R, K+1+2*N)*II*psi_t34*fi*pow(r, 2)
                fctp1 = gf*pow(r/R, K+2*N)*psi_t12*fi*pow(r, 2)
                fctp2 = ff*pow(r/R, K+2*N)*psi_t12*gi*pow(r, 2)

                intt1, intt2, intt3, intt4 = 0, 0, 0, 0
                intp1, intp2 = 0, 0
                dr = r[1] - r[0]
                for i in range(len(r)-1):
                    intt1 += (fctt1[i]+fctt1[i+1])/2 * dr
                    intt2 += (fctt2[i]+fctt2[i+1])/2 * dr
                    intt3 += (fctt3[i]+fctt3[i+1])/2 * dr
                    intt4 += (fctt4[i]+fctt4[i+1])/2 * dr
                    intp1 += (fctp1[i]+fctp1[i+1])/2 * dr
                    intp2 += (fctp2[i]+fctp2[i+1])/2 * dr
                
                ma += ME.MANKLs(N, K, L, s, k, m, n, rho, ki, kf, ni, nf, nui, nuf, Ji, r, R) * obd
                mt += -Constants.fT/R*np.sqrt(2)/np.sqrt(2*Ji+1)* (np.sqrt(K/(2*K+1))*(-np.sign(ki)*Utilities.GKLs(K, K, 1, kf, -ki)*intt1 + np.sign(kf)*Utilities.GKLs(K, K, 1, -kf, ki)*intt2) + (-Utilities.GKLs(K, K+1, 1, kf, ki)*intt3 + np.sign(ki)*np.sign(kf)*Utilities.GKLs(K, K+1, 1, -kf, -ki)*intt4)) * obd
                mp += Constants.fP/R*np.sqrt(2)/np.sqrt(2*Ji+1)*np.sqrt((K+1)/(2*K+1))*(-np.sign(ki)*Utilities.GKLs(K, K, 0, kf, -ki)*intp1 + np.sign(kf)*Utilities.GKLs(K, K, 0, -kf, ki)*intp2) * obd
               
        f.close()
    if beta_type == 'beta+' :
        fanklp1 = -(Constants.LAMBDA * float(ma) + float(mt) + float(mp))
    else :
        fanklp1 = -(Constants.LAMBDA * float(ma) - float(mt) + float(mp))
    return fanklp1
    
def FVNKLs(N, K, L, s, k, m, n, rho, chemin):
    if (s==0) :
        L = K
        FVNKLs = FVNKK0(N, K, L, s, k, m, n, rho, chemin)
    elif (s==1) :
        if (L==K) :
            FVNKLs = FVNKK1(N, K, L, s, k, m, n, rho, chemin)
        elif (L==K-1) :
            FVNKLs = FVNKLm1(N, K, L, s, k, m, n, rho, chemin)
        elif (L==K+1) :
            FVNKLs = FVNKLp1(N, K, L, s, k, m, n, rho, chemin)
    return FVNKLs

def FANKLs(N, K, L, s, k, m, n, rho, chemin):
    if (s==0) :
        L = K
        FANKLs = FANKK0(N, K, L, s, k, m, n, rho, chemin)
    elif (s==1) :
        if (L==K) :
            FANKLs = FANKK1(N, K, L, s, k, m, n, rho, chemin)
        elif (L==K-1) :
            FANKLs = FANKLm1(N, K, L, s, k, m, n, rho, chemin)
        elif (L==K+1) :
            FANKLs = FANKLp1(N, K, L, s, k, m, n, rho, chemin)
    return FANKLs

"""
chemin_210Bi = "/home/dumenil/Documents/THESE/DATA/OBTDs/test2.txt"
chemin_210Bi = "/home/dumenil/Documents/THESE/DATA/OBTDs/test.txt"

print("V^F(N=0, K=1, L=1, s=0, k=1, m=0, n=0, rho=0) =", FVNKLs(0, 1, 1, 0, 1, 0, 0, 0, chemin_210Bi))
print("V^F(N=0, K=1, L=1, s=0, k=1, m=1, n=1, rho=1) =", FVNKLs(0, 1, 1, 0, 1, 1, 1, 1, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=0, n=0, rho=0) =", FVNKLs(1, 1, 1, 0, 1, 0, 0, 0, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=1, n=1, rho=1) =", FVNKLs(1, 1, 1, 0, 1, 1, 1, 1, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=2, n=2, rho=1) =", FVNKLs(1, 1, 1, 0, 1, 2, 2, 1, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=2, n=2, rho=2) =", FVNKLs(1, 1, 1, 0, 1, 2, 2, 2, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=3, n=2, rho=2) =", FVNKLs(1, 1, 1, 0, 1, 3, 2, 2, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=3, n=3, rho=1) =", FVNKLs(1, 1, 1, 0, 1, 3, 3, 1, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=3, n=3, rho=2) =", FVNKLs(1, 1, 1, 0, 1, 3, 3, 2, chemin_210Bi))
print("V^F(N=1, K=1, L=1, s=0, k=1, m=3, n=3, rho=3) =", FVNKLs(1, 1, 1, 0, 1, 3, 3, 3, chemin_210Bi))
print("V^F(N=0, K=1, L=0, s=1, k=1, m=0, n=0, rho=0) =", FVNKLs(0, 1, 0, 1, 1, 0, 0, 0, chemin_210Bi))

print("A^F(N=0, K=1, L=1, s=1, k=1, m=0, n=0, rho=0) =", FANKLs(0, 1, 1, 1, 1, 0, 0, 0, chemin_210Bi))
print("A^F(N=0, K=1, L=1, s=1, k=1, m=1, n=1, rho=1) =", FANKLs(0, 1, 1, 1, 1, 1, 1, 1, chemin_210Bi))
print("A^F(N=1, K=1, L=1, s=1, k=1, m=0, n=0, rho=0) =", FANKLs(1, 1, 1, 1, 1, 0, 0, 0, chemin_210Bi))
print("A^F(N=1, K=1, L=1, s=1, k=1, m=1, n=1, rho=1) =", FANKLs(1, 1, 1, 1, 1, 1, 1, 1, chemin_210Bi))
"""

