#!/usr/bin/env python
import pysam
import sys
import numpy
import pandas as pd

if len(sys.argv) < 3 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
    print'''
    Created on Sep 19, 2013
    Updated by JEB on Oct 28, 2016

    Usage: Calculate_rpkm.py  accepted_hits_wt.bam accepted_hits_ko.bam output_prefix
    Note that the rpkm  is calculated here for the whole mRNA, not just the exons.

    2 additional functions were added on Nov 11, 2016:
    It is now possible to use your own gff3 or gtf file by adding the file name after the output_prefix argument
    e.g. Calculate_rpkm.py  accepted_hits_wt.bam accepted_hits_ko.bam output_prefix gtf_file

    It is also now possible to count reads in intergenic regions. Simply add "intergenic" as the last argument.
    e.g. Calculate_rpkm.py  accepted_hits_wt.bam accepted_hits_ko.bam output_prefix "intergenic"
    e.g. Calculate_rpkm.py  accepted_hits_wt.bam accepted_hits_ko.bam output_prefix gff3_file "intergenic"

    @author: chomer
    '''
    sys.exit()

'''populate annotation dictionary'''

if len(sys.argv) > 4:
    print sys.argv[4][-3:]
chrom_set = set()
if len(sys.argv) == 4 or sys.argv[4][-3:] not in ['ff3', 'gtf', 'gff']:
    gene_dict = {}
    fin = open("/home/jordan/GENOMES/CNA3_all_transcripts.gff3", "r")
    for line in fin:
        data = line.split("\t")
        if len(data) > 1 and data[2] == "mRNA":
            cnag = data[8].split(";")[0].split("=")[-1][0:10]
            chrom = data[0]
            chrom_set.add(chrom)
            start = data[3]
            end = data[4]
            strand = data[6]
            gene_dict[cnag] = [cnag, chrom, start, end, strand]

else:
    gene_dict = {}
    with open(sys.argv[4], 'r') as f:
        if sys.argv[4].endswith('gff3'):
            for line in f:
                data = line.split('\t')
                if len(data) > 1 and data[2] == 'mRNA':
                    cnag = data[8].split(";")[0].split("=")[-1][0:10]
                    chrom = data[0]
                    chrom_set.add(chrom)
                    start = data[3]
                    end = data[4]
                    strand = data[6]
                    gene_dict[cnag] = [cnag, chrom, start, end, strand]
        
        elif sys.argv[4].endswith('gtf'):
            for line in f:
                data = line.split('\t')
                if len(data) > 1 and data[2] == 'transcript':
                    cnag = data[8].split(';')[0].split('"')[1]
                    chrom = data[0]
                    chrom_set.add(chrom)
                    start = data[3]
                    end = data[4]
                    strand = data[6]
                    gene_dict[cnag] = [cnag, chrom, start, end, strand]
                    
print chrom_set
if 'intergenic' in sys.argv:
    print "Counting reads in intergenic regions"
    inter_dict = {}
    for chrom in list(chrom_set):
        chrom_genes = dict((k, gene_dict[k]) for k in gene_dict if gene_dict[k][1] == chrom)
        chr_gene_df = pd.DataFrame.from_dict(chrom_genes, orient='index')
        chr_gene_df.sort_values([0], inplace=True)
        sorted_genes = chr_gene_df.index.tolist()
        
        n = 0
        for n in range(len(sorted_genes)-1):
            gene = sorted_genes[n]
            next_gene = sorted_genes[n+1]
            gene_end = int(chr_gene_df[3][gene])
            next_start = int(chr_gene_df[2][next_gene])
            
            if next_start > gene_end:
                inter_name = gene+'_'+next_gene+'_intergenic'
                inter_dict[inter_name] = [inter_name, chrom, gene_end, next_start]

            else:
                print 'Overlapping transcripts:'
                print gene
                print next_gene
                
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

    iter_ko = fsam_ko.fetch(gene_data[1], int(gene_data[2]), int(gene_data[3]))
    for read in iter_ko:
        if gene_data[4] == "+":
            if read.is_reverse:
                ko_gene_count += 1
        elif gene_data[4] == "-":
            if not read.is_reverse:
                ko_gene_count += 1
                
    wt_gene_rpkm = float(wt_gene_count) * float(1000000) / (float(gene_length)*float(wt_reads))
    ko_gene_rpkm = float(ko_gene_count) * float(1000000) / (float(gene_length)*float(ko_reads))
    wt_rpkm_list.append(wt_gene_rpkm)
    wt_gene_list.append([gene, wt_gene_rpkm])
    ko_rpkm_list.append(ko_gene_rpkm)
    ko_gene_list.append([gene, ko_gene_rpkm])

    if int(ko_gene_count) == 0:
        ratio = sys.maxint
    else:
        ratio = float(wt_gene_rpkm) / float(ko_gene_rpkm)
    ratio2 = (float(wt_gene_count+1)*float(wt_reads)) / (float(ko_gene_count+1)*float(ko_reads))
    rpkm_dict[gene] = [wt_gene_rpkm, ko_gene_rpkm, float(ratio), wt_gene_count, ko_gene_count, ratio2]

    fsam_wt.close()
    fsam_ko.close()        

wt_inter_rpkms = []
wt_inter_list = []
ko_inter_rpkms = []
ko_inter_list = []

if 'intergenic' in sys.argv:
    for inter in inter_dict:
        fsam_wt = pysam.Samfile(sys.argv[1])
        fsam_ko = pysam.Samfile(sys.argv[2])
        wt_inter_count = 0
        ko_inter_count = 0
        inter_data = inter_dict[inter]
        inter_length = int(inter_data[3]) - int(inter_data[2])
        
        iter_wt = fsam_wt.fetch(inter_data[1], int(inter_data[2]), int(inter_data[3]))
        for read in iter_wt:
            wt_inter_count += 1
            
        iter_ko = fsam_ko.fetch(inter_data[1], int(inter_data[2]), int(inter_data[3]))
        for read in iter_ko:
            ko_inter_count += 1
            
        wt_inter_rpkm = float(wt_inter_count) * float(1000000) / (float(inter_length)*float(wt_reads))
        ko_inter_rpkm = float(ko_inter_count) * float(1000000) / (float(inter_length)*float(ko_reads))
        wt_inter_rpkms.append(wt_inter_rpkm)
        wt_inter_list.append([inter, wt_inter_rpkm])
        ko_inter_rpkms.append(ko_inter_rpkm)
        ko_inter_list.append([inter, ko_inter_rpkm])
    
    
        if int(ko_inter_count) == 0:
            ratio = sys.maxint
        else:
            ratio = float(wt_inter_rpkm) / float(ko_inter_rpkm)
        ratio2 = (float(wt_inter_count+1)*float(wt_reads)) / (float(ko_inter_count+1)*float(ko_reads))
        rpkm_dict[inter] = [wt_inter_rpkm, ko_inter_rpkm, float(ratio), wt_inter_count, ko_inter_count, ratio2]

        fsam_wt.close()
        fsam_ko.close()        
#print rpkm_dict
wt_median = numpy.median(wt_rpkm_list)
ko_median = numpy.median(ko_rpkm_list)
print wt_median
print ko_median

fout = open("{0}_rpkm_ratios.txt".format(sys.argv[3]), "w")
header_list = ["Gene", "WT RPKM", "KO RPKM", "WT/KO", "Wild Reads", "KO Reads", "Normalized WT/KO", "\n"]
#header_list = ["Gene", "WT RPKM", "KO RPKM", "WT/KO", "Wild Reads", "KO Reads", "\n"]
header = "\t".join(header_list)
fout.write(header)
ko_median_gene = ""
wt_median_gene = ""
ko_diff = sys.maxint
wt_diff = sys.maxint
for gene in rpkm_dict:
    data = rpkm_dict[gene]
    if wt_median - data[0] < wt_diff:
        wt_diff = wt_median - data[0]
        wt_median_gene = gene
    if ko_median - data[1] < ko_diff:
        ko_diff = ko_median - data[1]
        ko_median_gene = gene
    line_list = [gene, str(data[0]), str(data[1]), str(data[2]), str(data[3]), str(data[4]), str(data[5]), "\n"]
    #line_list = [gene, str(data[0]), str(data[1]), str(data[2]), str(data[3]), str(data[4]), "\n"]
    line = "\t".join(line_list)
    fout.write(line)

print wt_median_gene
print ko_median_gene

fout.close()
