'''Script for calling 5' and 3' splice sites and relative strengths of each sites 

Usage: python Call_splice_site_motifs.py

This will generate an output spreadsheet with relative splice-site strength for each transcript in the genome'''

import sys
import math
import numpy as np
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
import GeneUtility

def mirrored(listtree):
    return [(mirrored(x) if isinstance(x, list) else x) for x in reversed(listtree)]

#Populate gene dictionary and build genome
gene_list_as_dict = GeneUtility.GeneDict()
gene_list_as_dict.PopulateFromFile_new()
genome = GeneUtility.BuildGenome()
print genome.keys()

#First generate a consensus matrix for 5' and 3' splice site, where 1st row is A counts, second row is C, third row is T, fourth row is G.
pos_matrix_5prime = np.zeros([4,8])
pos_matrix_3prime = np.zeros([4,8])

for cnag in gene_list_as_dict:
    gene = gene_list_as_dict.FindGeneByCNAG(cnag)
    introns = gene.introns
    if gene.strand == "-":
        introns = mirrored(introns)
    for intron in introns:
        if gene.strand == "-":
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 7, intron[0] + 1, gene.strand, genome)
        else:
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 2, intron[0] + 6, gene.strand, genome)
        for a, base in enumerate(seq):
             if base == "A":
                 pos_matrix_5prime[0,a] = pos_matrix_5prime[0,a]+1
             if base == "C":
                 pos_matrix_5prime[1,a] = pos_matrix_5prime[1,a]+1
             if base == "T":
                 pos_matrix_5prime[2,a] = pos_matrix_5prime[2,a]+1
             if base == "G":
                 pos_matrix_5prime[3,a] = pos_matrix_5prime[3,a]+1

        if gene.strand == "-":
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[1] - 2, intron[1] + 6, gene.strand, genome)
            #print gene.strand
            #print seq
        else:
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[1] - 7, intron[1] + 1, gene.strand, genome)
            #print gene.strand
            #print seq
        for b, base in enumerate(seq):
            if base == "A":
                pos_matrix_3prime[0,b] = pos_matrix_3prime[0,b]+1
            if base == "C":
                pos_matrix_3prime[1,b] = pos_matrix_3prime[1,b]+1
            if base == "T":
                pos_matrix_3prime[2,b] = pos_matrix_3prime[2,b]+1
            if base == "G":
                pos_matrix_3prime[3,b] = pos_matrix_3prime[3,b]+1

a = 0
while a < 4:
    b = 0
    while b < 8:
        pos_matrix_5prime[a,b] = (pos_matrix_5prime[a,b])/36855.
        pos_matrix_3prime[a,b] = (pos_matrix_3prime[a,b])/36855.
        b += 1
    a += 1

print sum(pos_matrix_5prime)        
print pos_matrix_5prime
print sum(pos_matrix_3prime)
print pos_matrix_3prime

#Compare each splice site and output scores into a spreadsheet
fout = open("Splice_site_strengths.csv", "w")
header_list = ["Transcript","Intron", "5' Splice Site Score", "3' Splice Site Score", "Intron Length", "\n"]
header = "\t".join(header_list)
fout.write(header)

for cnag in gene_list_as_dict: 
    gene = gene_list_as_dict.FindGeneByCNAG(cnag)
    introns = gene.introns
    if gene.strand =="-":
        introns = mirrored(introns)
    for intron_num, intron in enumerate(introns):
        gene_matrix_5prime = np.zeros([4,8])
        if gene.strand == "-":
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 7, intron[0] + 1, gene.strand, genome)
        else:
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 2, intron[0] + 6, gene.strand, genome)
        for a, base in enumerate(seq):
             if base == "A":
                 gene_matrix_5prime[0,a] = gene_matrix_5prime[0,a]+1
             if base == "C":
                 gene_matrix_5prime[1,a] = gene_matrix_5prime[1,a]+1
             if base == "T":
                 gene_matrix_5prime[2,a] = gene_matrix_5prime[2,a]+1
             if base == "G":
                 gene_matrix_5prime[3,a] = gene_matrix_5prime[3,a]+1

        gene_matrix_3prime = np.zeros([4,8])
        if gene.strand == "-":
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 2, intron[0] + 6, gene.strand, genome)
            intron_seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[1], intron[0], gene.strand, genome)
        else:
            seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 7, intron[0] + 1, gene.strand, genome)
            intron_seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0], intron[1], gene.strand, genome)
        for b, base in enumerate(seq):
            if base == "A":
                gene_matrix_3prime[0,b] = gene_matrix_3prime[0,b]+1
            if base == "C":
                gene_matrix_3prime[1,b] = gene_matrix_3prime[1,b]+1
            if base == "T":
                gene_matrix_3prime[2,b] = gene_matrix_3prime[2,b]+1
            if base == "G":
                gene_matrix_3prime[3,b] = gene_matrix_3prime[3,b]+1

         #Calculate Scores (score of 0 is perfect, higher score is worse, 40 is highest possible)
        score_5prime = 0
        score_3prime = 0
        a = 0
        while a < 4:
            b = 0
            while b < 8:
                score_5prime += abs(pos_matrix_5prime[a,b] - (gene_matrix_5prime[a,b]))
                score_3prime += abs(pos_matrix_3prime[a,b] - (gene_matrix_3prime[a,b]))
                b += 1
            a += 1

        intron_length = len(intron_seq)

#        format_transcript_list = [cnag, "T0-", str(intron_num)]
#        format_transcript = "".join(format_transcript_list)
        lineout_list = [cnag+"T0", str(intron_num), str(score_5prime), str(score_3prime), str(intron_length), "\n"]
        lineout = "\t".join(lineout_list)
        fout.write(lineout)

fout.close()
