'''
Created on Oct 15, 2013

Usage: python Calculate_rpkm.py  accepted_hits_wt.bam accepted_hits_ko.bam output_prefix
Note that the rpkm  is calculated here for intergenic regions upstream of each gene. 
@author: jburke
'''

import pysam
import sys
import numpy
import GeneUtility

'''Creating a dictionary of gene objects, then a dictionary of intergenic regions referenced by CNAG'''
gene_dict = GeneUtility.GeneDict()
gene_dict.PopulateFromFile_new()
    
misc_dict = GeneUtility.MiscDict()
misc_dict.PopulateMiscFromFile()
misc_intergenic_dict = misc_dict.PopulateIntergenics()
    
intergenics_as_dict = GeneUtility.IntergenicDict()
intergenics_as_dict.PopulateIntergenicRegion(gene_dict)
intergenics_as_dict.FindTypes(gene_dict, misc_dict)

annotation_dict = GeneUtility.GeneAnnotationDict()
annotation_dict.PopulateGeneAnnotation_new()

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
print wt_reads

iter_ko = fsam_ko.fetch()
for read in iter_ko:
    if read.is_unmapped:
        pass
    else:
        ko_reads += 1
print ko_reads

fsam_wt.close()
fsam_ko.close()


'''compile read counts from wt and ko and put in dictionary'''
intergenic_dict = {}
wt_gene_list = []
ko_gene_list = []

for intergenic_CNAG in intergenics_as_dict:
    fsam_wt = pysam.Samfile(sys.argv[1])
    fsam_ko = pysam.Samfile(sys.argv[2])
    wt_gene_count = 0
    ko_gene_count = 0
    strand = gene_dict.FindGeneStrand(intergenic_CNAG)
    intergenic = gene_dict.FindIntergenicRegion(intergenic_CNAG, strand)
#    self.IntergenicsByCNAG[intergenic_CNAG] = intergenic
    chrom = "2.{0}".format(intergenic.chromosome)
    start = intergenic.start
    end = intergenic.end
    
    iter_wt = fsam_wt.fetch(chrom, int(start), int(end))
    for read in iter_wt:
        wt_gene_count += 1

    iter_ko = fsam_ko.fetch(chrom, int(start), int(end))
    for read in iter_ko:
        ko_gene_count += 1
                
    wt_gene_list.append(read)
    ko_gene_list.append(read)

    ratio2 = (float(wt_gene_count+1)*float(wt_reads)) / (float(ko_gene_count+1)*float(ko_reads)) 
    #To include reads with no counts in ko sample
    
    intergenic_dict[intergenic_CNAG] = [wt_gene_count, ko_gene_count, float(ratio2)]
    
#    print wt_gene_count
#    print ko_gene_count

    fsam_wt.close()
    fsam_ko.close()        

'''Write out tabulated file with read counts for wt and ko and the normalized ratio'''
fout = open("{0}_intergenic_read_ratios.txt".format(sys.argv[3]), "w")
header_list = ["Gene", "Wild Reads", "KO Reads", "Normalized Reads WT/KO", "\n"]
header = "\t".join(header_list)
fout.write(header)
ko_diff = sys.maxint
wt_diff = sys.maxint
for gene in intergenic_dict:
    data = intergenic_dict[gene]
    if int(data[0]) > 20 or int(data[1]) > 20:
        line_list = [gene, str(data[0]), str(data[1]), str(data[2]), "\n"]
        line = "\t".join(line_list)
        fout.write(line)

fout.close()
