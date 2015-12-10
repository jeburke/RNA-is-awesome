__author__ = 'jordanburke'

import sys

'''Usage: python LaSSO_alignment_db_cat.py <alignment_fasta_output> <LaSSO_intron_database>'''

fin1 = open(sys.argv[1], "r")
fin2 = open(sys.argv[2], "r")
fout = open("{0}_concatenated.fasta".format(sys.argv[1].split(".")[0]), "w")

intron_list = []
alignment_list = []

for j, line1 in enumerate(fin1):
    if j % 2 == 0:
        intron_list.append(line1.strip())
    elif j % 2 == 1:
        alignment_list.append(line1.strip())

alignment_di = dict(zip(intron_list, alignment_list))

intron_list_db = []
seq_list_db = []

for k, line2 in enumerate(fin2):
    if k % 2 == 0:
        intron_list_db.append(line2.strip())
    if k % 2 == 1:
        seq_list_db.append(line2.strip())

intron_db_di = dict(zip(intron_list_db, seq_list_db))

for intron, alignment in alignment_di.iteritems():
    fout.write(intron.strip()+"\n")
    fout.write(alignment.strip()+"\n")
    fout.write(intron_db_di[intron].strip()+"\n")

fin1.close()
fin2.close()
fout.close()