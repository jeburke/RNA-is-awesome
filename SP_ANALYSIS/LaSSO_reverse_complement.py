__author__ = 'jordanburke'

import sys

fin = open(sys.argv[1], "r")
fout = open("{0}_reversed.db".format(sys.argv[1].split(".")[0]), "w")

alt_map = {'ins':'0'}
complement = {'A':'T','C':'G','G':'C','T':'A'}

def reverse_complement(seq):
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases

for k, line in enumerate(fin):
    if k % 2 == 0:
        print line
        fout.write(line)
    elif k % 2 == 1:
        print line
        line = reverse_complement(line.strip())+"\n"
        print line
        fout.write(line)

fin.close()
fout.close()