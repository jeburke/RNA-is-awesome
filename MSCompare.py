'''Created on 20 May 2015 by Jordan E. Burke
This script can compare two lists of MS data and return genes present in both lists that have at least 10 counts in either replicate. Format input file as follows:
Column 1: Accession ID
Column 2: Counts Replicate 1
Column 3: Counts Replicate 2
Column 4: Coverage Replicate 1
Column 5: Coverage Replicate 2
Usage: python MSCompare.py <list1> <list2>
'''

import re
import sys
from collections import OrderedDict

def buildGeneList(file):
    genes = OrderedDict()
    with open(file, "r") as fin:
        for line in fin:
            if not line.startswith('C'):
                continue
            lineString = line.rstrip('\r\n')
            columns = re.split(r'\t+', lineString)
            genes[columns[0]] = (
                columns[1],  # Count A
                columns[2],  # Count B
                columns[3],  # Coverage A
                columns[4])  # Coverage B
    return genes

def pickGenes(genes):
    return OrderedDict(
        (g, (countA, countB, covA, covB))
        for g, (countA, countB, covA, covB) in genes.items()
        if max(countA, countB) > 10)

geneDict1 = buildGeneList(sys.argv[1])
geneDict2 = buildGeneList(sys.argv[2])

print len(geneDict1)
print len(geneDict2)

sGenes1 = pickGenes(geneDict1)
sGenes2 = pickGenes(geneDict2)

print len(sGenes1)
print len(sGenes2)

commonGenes = OrderedDict(
    (g, (count1A, count1B, sGenes2[g][0], sGenes2[g][1],
         cov1A, cov1B, sGenes2[g][2], sGenes2[g][3],
         (float(sGenes2[g][0]) + float(sGenes2[g][1])) / (float(count1A) + float(count1B))))
    for g, (count1A, count1B, cov1A, cov1B) in sGenes1.items()
    if g in sGenes2)
print len(commonGenes)

fmt_args = {
    'a_file': sys.argv[1].split('.')[0],
    'b_file': sys.argv[2].split('.')[0],
}

fout = open("{a_file}_{b_file}_comparison.txt".format(**fmt_args), "w")

header_list = ("Accession",
               "Counts {a_file}A", "Counts {a_file}B",
               "Counts {b_file}A", "Counts {b_file}B",
               "Coverage {a_file}A", "Coverage {a_file}B",
               "Coverage {b_file}A", "Coverage {b_file}B",
               "Ratio of average counts from {b_file}/{a_file}", "\n")
header = "\t".join(h.format(**fmt_args) for h in header_list)
fout.write(header)

for g, cols in commonGenes.items():
    line_list = list(g, *cols) + ["\n"]
    line = "\t".join(str(col) for col in line_list)
    fout.write(line)
