'''
Created on Jan 24, 2014

@author: jordan
'''

import sys

fin = open(sys.argv[1], "r")
fout = open("{0}_out_strip6.fastq".format(sys.argv[1].split(".")[0]), "w")

TDNA = "TGTCTAAGCGTCAATTTGTTTACACCACAATATATC"
T_min = 5

left = ""
headerline = ""

for k, line in enumerate(fin):
    if k % 4 == 1:
        left = 0
        for i in range(len(TDNA)):
            if line[i] != TDNA[i]:
                left = i
                break
        if left >= T_min:
            match = True
            fout.write(headerline)
            fout.write(line[left:])
    elif k%4 == 2 and match:
        fout.write(line)
    elif k%4 == 3 and match:
        fout.write(line[left:])
    elif k%4 == 0:
        match = False
        # cache for later, in case it *is* a match
        headerline = line


fin.close()
fout.close()
