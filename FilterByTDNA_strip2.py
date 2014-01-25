'''
Created on Jan 24, 2014

@author: jordan
'''

import sys

fin = open(sys.argv[1], "r")
fout = open("{0}_out_strip6.fastq".format(sys.argv[1].split(".")[0]), "w")

left = ""
headerline = ""
for k, line in enumerate(fin):
    n = 36
    TDNA = "TGTCTAAGCGTCAATTTGTTTACACCACAATATATC"
    if k % 4 == 1:
        left = 0
        while line[0:n] == TDNA:
            out = "yes"
            line = line[n:]
            left = n
            fout.write(headerline)
            fout.write(line)
            n = n - 1
            TDNA = TDNA[0:n]
    elif k%4 == 2 and out == "yes":
        fout.write(line)
    elif k%4 == 3 and out == "yes":
        fout.write(line[left:])
    elif k%4 == 0:
        out = "no"
        headerline = line


fin.close()
fout.close()
