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
    geneID = []
    coverageA = []
    coverageB = []
    countsA = []
    countsB = []
    fin = open(file, "r")
    for line in fin:
#        if line.startswith('P'):
        columns = re.split(r'\t+', line)
        geneID.append(columns[0])
        countsA.append(columns[1])
        countsB.append(columns[2])
        coverageA.append(columns[3])
        coverageB.append(columns[4])
    return geneID, coverageA, coverageB, countsA, countsB
    fin.close()

def pickGenes(geneList, countListA, countListB, coverageListA, coverageListB):
    n = 1
    selectGenes = []
    selectCountsA = []
    selectCountsB = []
    selectCoverageA = []
    selectCoverageB = []
    while n < len(geneList):
        if int(countListA[n]) > 10 or int(countListB[n]) > 10:
            selectGenes.append(geneList[n])
            selectCountsA.append(countList[n])
            selectCountsB.append(countList[n])
            selectCoverageA.append(coverageList[n])
            selectCoverageB.append(coverageList[n])
        n += 1
    return selectGenes, selectCountsA, selectCountsB, selectCoverageA, selectCoverageB

geneList1, coverageList1a, coverageList1b, countsList1a, countsList1b = buildGeneList(sys.argv[1])
geneList2, coverageList2a, coverageList2b, countsList2a, countsList2b = buildGeneList(sys.argv[2])

print len(geneList1)
print len(geneList2)

sGenes1, sCounts1a, sCounts1b, sCoverage1a, sCoverage1b = pickGenes(geneList1, countsList1a, countsList1b, coverageList1a, coverageList1b)
sGenes2, sCounts2a, sCounts2b, sCoverage2a, sCoverage2b = pickGenes(geneList2, countsList2a, countsList2b, coverageList2a, coverageList2b)

print len(sGenes1)
print len(sGenes2)

finalGenes1 = []
finalCounts1a = []
finalCounts1b = []
finalCounts2a = []
finalCounts2b = []
finalCoverage1 = []
for i in range(len(sGenes1)):
    if sGenes1[i] in sGenes2:
        finalGenes1.append(sGenes1[i])
        finalCounts1a.append(sCounts1a[i])
        finalCounts1b.append(sCounts1b[i])
        finalCounts2a.append(sCounts2a[i])
        finalCounts2b.append(sCounts2b[i])
        finalCoverage1a.append(sCoverage1a[i])
        finalCoverage1b.append(sCoverage1b[i])
        finalCoverage2a.append(sCoverage2a[i])
        finalCoverage2b.append(sCoverage2b[i])

#finalList2 = []
#for i in range(len(geneList2)):
#    if geneList2[i] not in geneList1:
#        finalList2.append(geneList2[i])

a=sys.argv[1].split(".")
b=sys.argv[2].split(".")

fout = open("{0}_{1}_comparison.txt".format(a[0], b[0]), "w")

#fout.write("Accession","\t","Counts 1A","\t","Counts 1B","\t","Counts2A","\t","Counts 2B","\t","Coverage 1A","\t","Coverage 1B","\t","Coverage 2A","\t","Coverage 2B","\n")
for i in range(len(finalGenes1)):
    line_list = [str(finalGenes1[i]), str(finalCounts1a[i]), str(finalCounts1b[i]), str(finalCounts2a[i]), str(finalCounts2b[i]), str(finalCoverage1a[i]), str(finalCoverage1b[i]), str(finalCoverage2a[i]), str(finalCoverage2b[i]), "\n"]
    line = "\t".join(line_list)
    fout.write(line)

