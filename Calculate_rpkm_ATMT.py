'''
Created on Sep 19, 2013

Usage: python Calculate_rpkm.py  accepted_hits_wt.bam accepted_hits_ko.bam output_prefix
Note that the rpkm  is calculated here for the whole mRNA, not just the exons. 
@author: chomer
'''

import pysam
import sys
import numpy


'''populate annotation dictionary'''
gene_dict = {}
fin = open("/home/chomer/Sequencing/CNA2_Combined.gff3", "r")
for line in fin:
    data = line.split("\t")
    if len(data) > 1 and data[2] == "mRNA":
        cnag = data[8].split(";")[0].split("=")[-1][0:10]
        chrom = "2.{0}".format(data[0].split(".")[-1])
        start = data[3]
        end = data[4]
        strand = data[6]
        gene_dict[cnag] = [cnag, chrom, start, end, strand]


'''determine the number of reads in each dataset'''
fsam_wt = pysam.Samfile(sys.argv[1])
fsam_ko = pysam.Samfile(sys.argv[2])
wt_reads = 0
ko_reads = 0
iter_wt = fsam_wt.fetch()
for read in iter_wt:
    if read.is_unmapped:
        pass
    else:
        wt_reads += 1
iter_ko = fsam_ko.fetch()
for read in iter_ko:
    if read.is_unmapped:
        pass
    else:
        ko_reads += 1
fsam_wt.close()
fsam_ko.close()
print wt_reads
print ko_reads


'''compile read counts from wt and ko and put in dictionary'''
rpkm_dict = {}
wt_rpkm_list = []
ko_rpkm_list = []
wt_gene_list = []
ko_gene_list = []

for gene in gene_dict:
    fsam_wt = pysam.Samfile(sys.argv[1])
    fsam_ko = pysam.Samfile(sys.argv[2])
    wt_gene_count = 0
    ko_gene_count = 0
    gene_data = gene_dict[gene]
    gene_length = int(gene_data[3]) - int(gene_data[2])

    iter_wt = fsam_wt.fetch(gene_data[1], int(gene_data[2]), int(gene_data[3]))
    for read in iter_wt:
        if gene_data[4] == "+":
            if read.is_reverse:
                wt_gene_count += 1
        elif gene_data[4] == "-":
            if not read.is_reverse:
                wt_gene_count += 1
        elif gene_data[4] == "-":
            if read.is_reverse:
                wt_gene_count += 1
        elif gene_data[4] == "+":
            if not read.is_reverse:
                wt_gene_count += 1
#        wt_gene_count += 1

    iter_ko = fsam_ko.fetch(gene_data[1], int(gene_data[2]), int(gene_data[3]))
    for read in iter_ko:
        if gene_data[4] == "+":
            if read.is_reverse:
                ko_gene_count += 1
        elif gene_data[4] == "-":
            if not read.is_reverse:
                ko_gene_count += 1
        elif gene_data[4] == "-":
            if read.is_reverse:
                ko_gene_count += 1
        elif gene_data[4] == "+":
            if not read.is_reverse:
                ko_gene_count += 1

    wt_gene_rpkm = float(wt_gene_count) * float(1000000) / (float(gene_length)*float(wt_reads))
    ko_gene_rpkm = float(ko_gene_count) * float(1000000) / (float(gene_length)*float(wt_reads))
    wt_rpkm_list.append(wt_gene_rpkm)
    wt_gene_list.append([gene, wt_gene_rpkm])
    ko_rpkm_list.append(ko_gene_rpkm)
    ko_gene_list.append([gene, ko_gene_rpkm])

    if int(ko_gene_count) == 0:
        ratio = sys.maxint
    else:
        ratio = float(wt_gene_rpkm) / float(ko_gene_rpkm)
#    if int(ko_gene_count) == 0: #To exclude reads with no counts in ko sample
#        ratio2 = "NA"
#    else:
#        ratio2 = (float(wt_gene_count)*float(wt_reads)) / (float(ko_gene_count)*float(ko_reads))
    ratio2 = (float(wt_gene_count+1)*float(ko_reads)) / (float(ko_gene_count+1)*float(wt_reads)) #To include reads with no counts in ko sample
    rpkm_dict[gene] = [wt_gene_rpkm, ko_gene_rpkm, float(ratio), wt_gene_count, ko_gene_count, float(ratio2)]

    fsam_wt.close()
    fsam_ko.close()        

#wt_median = numpy.median(wt_rpkm_list)
#ko_median = numpy.median(ko_rpkm_list)
#print wt_median
#print ko_median
fout = open("{0}_rpkm_ratios.txt".format(sys.argv[3]), "w")
header_list = ["Gene", "WT RPKM", "KO RPKM", "RPKM WT/KO", "Wild Reads", "KO Reads", "Normalized Reads WT/KO", "\n"]
header = "\t".join(header_list)
fout.write(header)
#ko_median_gene = ""
#wt_median_gene = ""
#ko_diff = sys.maxint
#wt_diff = sys.maxint
for gene in rpkm_dict:
    data = rpkm_dict[gene]
#    if wt_median - data[0] < wt_diff:
#        wt_diff = wt_median - data[0]
#        wt_median_gene = gene
#    if ko_median - data[1] < ko_diff:
#        ko_diff = ko_median - data[1]
#        ko_median_gene = gene
#    line_list = [gene, str(data[0]), str(data[1]), str(data[2]), str(data[3]), str(data[4]), str(data[2]/wt_median*ko_median), "\n"]
    if int(data[3]) > 20 or int(data[4]) > 20:
        line_list = [gene, str(data[0]), str(data[1]), str(data[2]), str(data[3]), str(data[4]), str(data[5]), "\n"]
        line = "\t".join(line_list)
        fout.write(line)

#print wt_median_gene
#print ko_median_gene

fout.close()
