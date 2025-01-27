#################################################
# Read_OBTD.py file
# Read the OBTD file in kshell and translate this file to the calculations
# python3 Read_OBTD.py fileName
#################################################
import sys

fileName = "OBTD_R1_Bi210_Po210_khpe.dat" #sys.argv[1]

"""
betatype = sys.argv[2]
A = sys.argv[3]
Z = sys.argv[4]
Ji = sys.argv[5]
nucleus = str(sys.argv[7])
"""
"""
def k_identification(n, l, j):
    j = j/2
    if l == 0 : # s
        return -1, n, l
    elif l == 1 : # p
        if j == 1/2 :
            return 1, n, l
        elif j == 3/2 :
            return -2, n, l
    elif l == 2 : # d
        if j == 3/2 :
            return 2, n, l
        elif j == 5/2 :
            return -3, n, l
    elif l == 3 : # f
        if j == 5/2 :
            return 3, n, l
        elif j == 7/2 :
            return -4, n, l
    elif l == 4 : # g
        if j == 7/2 :
            return , n, l
        elif j == 9/2 : 
            return , n, l
    elif l == 5 : # h
        if j == 9/2 :
            return , n, l
        elif j == 11/2 :
            return , n, l 
    elif l == 6 : # i
        if j == 11/2 :
            return , n, l
        elif j == 13/2 :
            return , n, l
    elif l == 7 : # j
        if j == 13/2 :
            return , n, l
        elif j == 15/2 :
            return , n, l
"""            
def find_rank(file): # equivalent at delta J
    with open(file, newline='\n') as fichier_txt: 
        f = open(file, 'r')
        for lines in fichier_txt:
            data = lines.split()
            if len(data) == 5 :
                if data[1] == 'rank' :
                    rank = int(data[2])
                    print(rank) 
        f.close()
        return rank
                 
def data_transition(file) : 
    data_transition = []
    with open(file, newline='\n') as fichier_txt: 
        f = open(file, 'r')
        for lines in fichier_txt:
            data = lines.split()
            
            if (len(data) == 5) :
                if data[0] != '#' :
                    
                    data_transition.append(data) 
        
        del(data_transition[0])
        del(data_transition[0])
        del(data_transition[0])            
        print(data_transition)
        f.close()
        return data_transition

def data_index(file) :
    data_idx = []
    with open(file, newline='\n') as fichier_txt: 
        f = open(file, 'r')
        for lines in fichier_txt:
            data = lines.split()
            if len(data) == 6 :
                del(data[0])
                data_idx.append(data)
                #print(data)
        print(data_idx)
        f.close()

betatype = 'beta-'        
A = '210'
Z = '83'
dj = '1'
nucleus = 'Bi'
def translate(file) :
    data_transi = data_transition(file)
    data_idx = data_index(file)
    with open('OBTD_'+A+nucleus, 'w') as fichier_txt: 
        fichier_txt.write('#\t'+betatype+'\t'+A+'\t'+Z+'\t'+dj+'\t'+nucleus+'\t'+'#\n')
        for i in range(len(data_transi)):
            idx_i = data_transi[i][0]
            idx_j = data_transi[i][0]
            print(idx_i)
        
translate(fileName)        
#data_transition(fileName)
#find_rank(fileName)
#data_index(fileName)

