'''Created on 4 June 2014 by Jordan E. Burke
This script can compare any two lists and return items in one list but not the other. Listed one after the other in the output file.
Usage: python CompareLists.py <list1> <list2>
'''

import sys
from collections import OrderedDict

def buildList(file):
    geneID = []
    fin = open(file, "r")
    for line in fin:
        geneID.append(line.strip())
    return geneID
    fin.close()


geneList1 = buildList(sys.argv[1])
geneList2 = buildList(sys.argv[2])

finalList1 = []
for i in range(len(geneList1)):
    if geneList1[i] not in geneList2:
        finalList1.append(geneList1[i])

finalList2 = []
for i in range(len(geneList2)):
    if geneList2[i] not in geneList1:
        finalList2.append(geneList2[i])

a=sys.argv[1].split("_")
b=sys.argv[2].split("_")

fout = open("{0}_{1}_compared_list.txt".format(a[0], b[0]), "w")

fout.write("Only in {0} \n".format(a[0], b[0]))
for i in range(len(finalList1)):
    line_list = [str(finalList1[i]), "\n"]
    line = "\t".join(line_list)
    fout.write(line)

fout.write("\n Only in {1} \n".format(b[0], a[0]))
for i in range(len(finalList2)):
    line_list = [str(finalList2[i]), "\n"]
    line = "\t".join(line_list)
    fout.write(line)

fout.close()
