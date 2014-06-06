'''Created on 4 June 2014 by Jordan E. Burke

This script identifies genes or proteins in a blastn, blastp or blastx search 
that have no hits and returns a simple lest of gene identifiers.

Usage: python sortBlastn.py <BlastOutput> <FilePrefix> '''

import sys

geneID = []
output = []

def findHits(file):
    fin = open(file, "r")
    for line in fin:
        if line.startswith("Query="):
            geneID.append(line[7:18]) 
        elif line.startswith("Sequences") or line.startswith("*****"):
            output.append(line.strip())
    fin.close()

findHits(sys.argv[1])

print len(geneID)
print len(output)

fout = open("{0}_blastn_hits.txt".format(sys.argv[2]), "w")

#for i in range(len(geneID)):
#    line_list = [str(geneID[i]), str(output[i]), "\n"]
#    line = "\t".join(line_list)
#    fout.write(line)

for i in range(len(geneID)):
    if output[i] == "***** No hits found *****":
        line_list = [str(geneID[i]), "\n"]
        line = "\t".join(line_list)
        fout.write(line)

fout.close()
