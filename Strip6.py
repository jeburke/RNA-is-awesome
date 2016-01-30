#!/usr/bin/env python

''' Created on Aug 31, 2015
@author: jordan
'''

import sys

fin = open(sys.argv[1], "r")
if len(sys.argv) == 2:
    fout = open("{0}_debarcoded.fasta".format(sys.argv[1].split(".")[0]), "w")
elif len(sys.argv) == 3:
    fout = open("{0}_debarcoded.fasta".format(sys.argv[2]), "w")
else:
    print "Please provide input file"

left = ""
headerline = ""
for k, line in enumerate(fin):
    n = 6
    if k % 2 == 0:
#        print line
        fout.write(line)
    elif k % 2 == 1:
#        print line
        left = n
        line = line.rstrip()
        line = line[left:]+"\n"
        fout.write(line)

fin.close()
fout.close()
