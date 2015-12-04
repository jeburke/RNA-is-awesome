''' Created on Aug 31, 2015
@author: jordan
'''

import sys

fin = open(sys.argv[1], "r")
fout = open("{0}_LaSSO.fasta".format(sys.argv[1].split(".")[0]), "w")

left = ""
headerline = ""
for k, line in enumerate(fin):
    if k % 2 == 1:
        fout.write(line)

fin.close()
fout.close()
