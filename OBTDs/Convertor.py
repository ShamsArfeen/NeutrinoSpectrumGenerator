Sym="Bi"
Z=83
A=214
Ji=1.0
obd_files = ["OBTD_R1_Bi214_khpe_j2n_Po214_khpe_j0p.dat"] # all obd files to be coverted into a single input format

def convert_obd(obd_files,Z,A,Sym,Ji,file="OBTDS.txt"):

    f = open(file,'w+')
    f.write("# beta-\t"+str(A)+'\t'+str(Z)+'\t'+str(Ji)+'\t'+Sym+' #')

    for obd_file in obd_files:
        f2 = open(obd_file,"r")
        lines = f2.readlines()
        entries = []
        for i in range(2,len(lines)):
            entry = lines[i].split()
            if len(entry) > 0:
                entries.append((int(entry[2])+1, int(entry[3]), int(entry[4]), int(entry[5])))
            else:
                iline = i
                break
                
        rank = int(lines[iline+1].split()[2])
        parity = int(lines[iline+1].split()[4])
        no_of_obd = int(lines[iline+4].split()[4])

        for k in range(iline+11,iline+11+no_of_obd):
            f.write('\n')

            i = int(lines[k].split()[0])
            j = int(lines[k].split()[1])
            obd = float(lines[k].split()[2])

            ni = entries[i-1][0]
            nf = entries[j-1][0]
            
            ji = entries[i-1][2]/2
            jf = entries[j-1][2]/2
            
            li = entries[i-1][1]
            lf = entries[j-1][1]
            
            if ji==li-0.5:
                kappa_i = li
            elif ji==li+0.5:
                kappa_i=-(li+1)
            # ki = abs(kappa_i)

            if jf==lf-0.5:
                kappa_f= lf
            elif jf==lf+0.5:
                kappa_f=-(lf+1)
            # kf = abs(kappa_f)

            f.write(str(rank))
            f.write('\t')
            f.write(str(kappa_i))
            f.write('\t')
            f.write(str(kappa_f))
            f.write('\t')
            f.write(str(ni))
            f.write('\t')
            f.write(str(nf))
            f.write('\t')
            f.write("{:.6f}".format(obd))
    f.close()


convert_obd(obd_files,Z,A,Sym,Ji)

