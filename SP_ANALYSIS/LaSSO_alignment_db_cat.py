__author__ = 'jordanburke'

import sys

'''Usage: python LaSSO_alignment_db_cat.py <alignment_fasta_output> <LaSSO_intron_database>'''

fin1 = open(sys.argv[1], "r")
fin2 = open(sys.argv[2], "r")
fout1 = open("{0}_concatenated.fasta".format(sys.argv[1].split(".")[0]), "w")
fout2 = open("{0}_branch_points.txt".format(sys.argv[1].split(".")[0], "w")

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

header_list = ["Transcript","\t","Intron","\t","BP Position","\n"]
header = "".join(header_list)
fout2.write(header)

for intron, alignment in alignment_di.iteritems():
    if intron.split("_")[2] != "int":
        if int(intron.split("@")[1]) > 6:
            fout1.write(intron.strip()+"\n")
            fout1.write("Alignment:"+alignment.strip()+"\n")
            fout1.write("Reference:"+intron_db_di[intron].strip()+"\n"+"\n")

            line_list = [str(intron.split(";")[0]).strip(">"),"\t",str(intron.split(":")[2]),"\t",str(intron.split("@")[1]),"\n"]
            line = "".join(line_list)
            fout2.write(line_list)

fin1.close()
fin2.close()
fout1.close()
fout2.close())
