import radioactivedecay as rd
import requests
import csv

def queryDecayChain(nuc, population=1.0, decay={}):
    decay.setdefault(nuc, dict())
    # print(nuc)
    try:
      rdnuc = rd.Nuclide(nuc)
    except:
      return decay

    for d_mode, d_ratio, d_nuc in zip(rdnuc.decay_modes(), rdnuc.branching_fractions(), rdnuc.progeny()):
        # print(nuc, d_mode, d_ratio, d_nuc)
        if d_mode == 'β-':

            sym_A = nuc.split('-')
            if nuc[-1].isalpha():
                URL = "https://www-nds.iaea.org/relnsd/v1/data?fields=decay_rads&nuclides=" + sym_A[1][:-1] + sym_A[0].lower() + "&rad_types=bm"
                state = chr(ord(nuc[-1]) - ord('m') + ord('X'))
            else:
                URL = "https://www-nds.iaea.org/relnsd/v1/data?fields=decay_rads&nuclides=" + sym_A[1] + sym_A[0].lower() + "&rad_types=bm"
                state = ''

            r = requests.get(url = URL)
            data = r.text
            reader = csv.reader(data.strip().splitlines())

            sym_D = d_nuc.split('-')
            if d_nuc[-1].isalpha():
                d_URL = "https://www-nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=" + sym_D[1][:-1] + sym_D[0].lower()
                d_state = chr(ord(d_nuc[-1]) - ord('m') + ord('X'))
            else:
                d_URL = "https://www-nds.iaea.org/relnsd/v1/data?fields=levels&nuclides=" + sym_D[1] + sym_D[0].lower()
                d_state = ''
            r_d = requests.get(url = d_URL)
            data_d = r_d.text
            reader_d = csv.reader(data_d.strip().splitlines())

            max_energy = []
            intensity_beta = []
            p_energy = []
            transition = []
            p_j = []
            d_j = []
            for row in reader:
                if len(row) > 17 and row[0] != 'mean_energy' and row[2] != '' and row[6] != '' and row[17] != '' and row[16]==state:

                    max_energy.append(float(row[6]))
                    intensity_beta.append(float(row[2]))
                    p_energy.append(float(row[17]))
                    transition.append(row[10])
                    p_j.append(row[19])

                    daughter_energy_shift = ''
                    if row[5][-1].isalpha():
                        daughter_level_energy = float(row[5][:-2])
                        daughter_energy_shift = row[5][-1]
                    else:
                        daughter_level_energy = float(row[5])
                    diff = float('inf')
                    d_j_val = 'unknown'
                    reader_d = csv.reader(data_d.strip().splitlines())
                    for row_d in reader_d:
                        if  len(row_d) > 8 and row_d[5] != 'energy':
                            daughter_level_energy_d = float(row_d[5])
                            if abs(daughter_level_energy - daughter_level_energy_d) < diff and row_d[4] == daughter_energy_shift:
                                diff = abs(daughter_level_energy - daughter_level_energy_d)
                                d_j_val = row_d[8]
                    d_j.append(d_j_val)

            if len(p_energy) == 0:
                continue

            for i in range(len(max_energy)):
                decay[nuc].setdefault(max_energy[i], {})
                decay[nuc][max_energy[i]].setdefault('Contribution', 0)
                decay[nuc][max_energy[i]]['Contribution'] += population * intensity_beta[i] / 100
                decay[nuc][max_energy[i]]['Type'] = transition[i]
                decay[nuc][max_energy[i]]['Ji'] = p_j[i]
                decay[nuc][max_energy[i]]['Jf'] = d_j[i]


        queryDecayChain(d_nuc, population * d_ratio, decay)
    return decay

  