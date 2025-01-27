import numpy as np

def getJiJf(j_ini, j_fin):

    if len(j_ini) == 0 or len(j_fin) == 0:
        return None,None 

    if '+' in j_ini:
        par_i = '+'
    elif '-' in j_ini:
        par_i = '-'

    if '+' in j_fin:
        par_f = '+'
    elif '-' in j_fin:
        par_f = '-'

    j_ini = j_ini.split()[0].replace('(','').replace(')','').replace('+','').replace('-','')
    j_fin = j_fin.split()[0].replace('(','').replace(')','').replace('+','').replace('-','')

    if '/' in j_ini and '/' in j_fin:
        j_ini = float(j_ini.split('/')[0]) / float(j_ini.split('/')[1])
        j_fin = float(j_fin.split('/')[0]) / float(j_fin.split('/')[1])
    elif '/' not in j_ini and '/' not in j_fin:
        j_ini = float(j_ini)
        j_fin = float(j_fin)
    else:
        print('Bad values:', j_ini, j_fin)
        return None,None

    return j_ini, j_fin