__author__ = 'jordanburke'

import sys

'''Usage: python LaSSO_alignment_db_cat.py <alignment_fasta_output> <LaSSO_intron_database>'''

fin1 = open(sys.argv[1], "r")
fin2 = open(sys.argv[2], "r")
fout1 = open("{0}_concatenated.fasta".format(sys.argv[1].split(".")[0]), "w")
fout2 = open("{0}_branch_points.txt".format(sys.argv[1].split(".")[0]), "w")
fout3 = open("{0}_db_bps.txt".format(sys.argv[1].split(".")[0]), "w")

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

header_list = ["Intron","\t","BP Position","\n"]
header = "".join(header_list)
fout2.write(header)

position_list = []
for intron, alignment in alignment_di.iteritems():
    if intron.split("_")[2] != "int":
        if int(intron.split("@")[1]) > 6:
            fout1.write(intron.strip()+"\n")
            fout1.write("Alignment:"+alignment.strip()+"\n")
            fout1.write("Reference:"+intron_db_di[intron].strip()+"\n"+"\n")
            position_list.append(intron)

position_list_sorted = position_list.sort()

for position in position_list:
        line_list = [str(position.split(";")[0]).strip(">"),"-",str(str(position.split(":")[2]).split("-")[0]),"\t",str(position.split("@")[1]),"\n"]
        line = "".join(line_list)
        fout2.write(line)

fout3.write("Intron"+"\n")
for intron, reference in intron_db_di.iteritems():
    if intron.split("_")[2] != "int":
        if int(intron.split("@")[1]) > 6:
	    line_list = [str(intron.split(";")[0]).strip(">"),"-",str(str(position.split(":")[2]).split("-")[0]),"\n"]
            line = "".join(line_list)
            fout3.write(line)

fin1.close()
fin2.close()
fout1.close()
fout2.close()
