__author__ = 'jordanburke'

import pandas
import math
import random
import re
import numpy
import peakutils
import collections
from peakutils.plot import plot as pplot
from matplotlib import pyplot

##################################################################################
## Determines splice site locations from gff3 file. Needs to have "chr" format  ##
##################################################################################

def list_splice_sites(gff3_file, chromosome="All", gene_list=None):
    fin = open(gff3_file,"r")
    splice_site_dict = {}
    n = 1
    for line in fin:
        columns = re.split(r'\t+', line.strip())
        if len(columns) > 1:
            if columns[2] == "mRNA":
                CNAG = columns[8].strip()
                CNAG = CNAG[3:15]
                splice_site_dict[CNAG] = [[],[],columns[1]]
            elif columns[2] == "exon":
                CNAG = columns[8].strip()
                CNAG = CNAG[-12:]
                if gene_list is None:
                    if columns[6] == "+":
                        splice_site_dict[CNAG][0].append(int(columns[4])-1)
                        splice_site_dict[CNAG][1].append(int(columns[3])-2)
                    if columns[6] == "-":
                        splice_site_dict[CNAG][0].append(int(columns[3])-1)
                        splice_site_dict[CNAG][1].append(int(columns[4]))
                else:
                    if CNAG in gene_list:
                        if columns[6] == "+":
                            splice_site_dict[CNAG][0].append(int(columns[4])-1)
                            splice_site_dict[CNAG][1].append(int(columns[3])-2)
                        if columns[6] == "-":
                            splice_site_dict[CNAG][0].append(int(columns[3])-1)
                            splice_site_dict[CNAG][1].append(int(columns[4]))
    
    #Trim to entries in gene list
    if gene_list is not None:
        splice_site_dict = dict([(CNAG, splice_site_dict[CNAG]) for CNAG in gene_list])
    
    #Trim to just one chromosome
    if chromosome != "All":
        splice_site_dict = dict([(CNAG, coords) for CNAG, coords in splice_site_dict.iteritems() if coords[2] == chromosome])
    
    print len(splice_site_dict)
    return splice_site_dict
    
######################################################################################################
## Build a transcript dictionary with all information from a gff3 file                              ##
######################################################################################################

def build_transcript_dict(gff3_file):
    with open(gff3_file,"r") as gff3:
        transcript_dict = {}
        for line in gff3:
            columns = re.split(r'\t+', line)
            if len(columns) > 1:
                if columns[2] == "mRNA":
                    CNAG = columns[8]
                    CNAG = CNAG[3:15]
                    #Transcript dictionary: keys are CNAG, values are [start, end, strand, chromosome]
                    transcript_dict[CNAG] = [int(columns[3]), int(columns[4]), columns[6], columns[0]]
               
    transcript_dict = collections.OrderedDict(sorted(transcript_dict.items()))
    return transcript_dict
    
######################################################################################################
## Pick peaks using PeakUtils package. Can be further refined or substituted with another package   ##
######################################################################################################

def peaks_by_gene(gff3_file, bedgraph_file, chromosome="All", gene_list=None):
    #Create dictionary of known splice sites by transcript
    splice_site_dict = list_splice_sites(gff3_file, chromosome, gene_list=gene_list)
    
    #Build transcript dictionary
    transcript_dict = build_transcript_dict(gff3_file)
    
    #Initiate output files
    with open("{0}_{1}_indexes.txt".format(bedgraph_file.split("/")[-1].split(".")[0], chromosome),"w") as f_handle:
        line = "transcript\t"+"chromosome"+"peak type"+" coordinate\t"+"peak size\n"
        f_handle.write(line)
    with open("{0}_{1}_detection_rate.txt".format(bedgraph_file.split("/")[-1].split(".")[0], chromosome),"w") as f_handle:
        line = "transcript\t 5' splice site hits\t 5' splice sites \t 3' splice site hits \t 3' splice sites \t new peaks \n"
        f_handle.write(line)
    
    if gene_list is not None:
        transcript_dict = dict([(CNAG, transcript_dict[CNAG]) for CNAG in gene_list])
    
    if chromosome != "All":
        transcript_dict = dict([(CNAG, coords) for CNAG, coords in splice_site_dict.iteritems() if coords[3] == chromosome])
    
    fivep_hits = []
    threep_hits = []
    new_peaks = []
    counter = 0
    fivep_ss_total = 0
    threep_ss_total = 0
    
    #Read bedgraph and sort by transcript into dictionary                
    bedgraph_dict = {}
    for CNAG, value in transcript_dict.iteritems():
        bedgraph_dict[CNAG] = [[],[]]

    with open(bedgraph_file, "r") as bedgraph:
        for line in bedgraph:
            columns = re.split(r'\t', line)
            for CNAG, coords in transcript_dict.iteritems():
                #Empty dictionary for bedgraph. Values will be [list of genomic positions][reads starting at that position]
                if int(columns[1]) > coords[0] and int(columns[1]) < coords[1]:
                    bedgraph_dict[CNAG][0].append(int(columns[1]))
                    bedgraph_dict[CNAG][1].append(int(columns[3]))
    bedgraph_dict = collections.OrderedDict(sorted(bedgraph_dict.items()))
    
    #Create numpy arrays to pick peaks from bedgraph dictionary
    for CNAG, values in bedgraph_dict.iteritems():
        bygene_x = numpy.array(values[0])
        bygene_y = numpy.array(values[1])
        
        if numpy.sum(bygene_y) > 100:
            base = peakutils.baseline(bygene_y, 2)
            indexes = peakutils.indexes(bygene_y-base, thres=0.05, min_dist=5)
            
            if gene_list is not None:
                print CNAG
                pyplot.figure(figsize=(15,6))
                pplot(bygene_x, bygene_y-base, indexes)
            
            #Filter out peaks with less than 5 reads
            a = 0
            bygene_x_filt = []
            bygene_y_filt = []
            while a < len(bygene_x[indexes]):
                if bygene_y[indexes][a] > 5:
                    bygene_x_filt.append(bygene_x[indexes][a])
                    bygene_y_filt.append(bygene_y[indexes][a])
                a += 1
                
                                                          
            #Compare peaks to known splice sites from splice site dictionary                                             
            CNAG_fivep_hits = []
            CNAG_fivep_heights = []
            CNAG_threep_hits = []
            CNAG_threep_heights = []
            CNAG_new_peaks = []
            CNAG_new_peaks_heights = []
                                                          
            for dict_CNAG, sites in splice_site_dict.iteritems():
                if dict_CNAG == CNAG:
                    fivep_ss_total += len(sites[0])
                    threep_ss_total += len(sites[1])
                    i = 0
                    while i < len(bygene_x_filt):
                        if bygene_x_filt[i] in sites[0]:
                            CNAG_fivep_hits.append(bygene_x_filt[i])
                            CNAG_fivep_heights.append(bygene_y_filt[i])
                            i += 1
                        elif bygene_x_filt[i] in sites[1]:
                            CNAG_threep_hits.append(bygene_x_filt[i])
                            CNAG_threep_heights.append(bygene_y_filt[i])
                            i += 1
                        else:
                            CNAG_new_peaks.append(bygene_x_filt[i])
                            CNAG_new_peaks_heights.append(bygene_y_filt[i])
                            i += 1
                                                          
            #Add number for this transcript to existing totals
            fivep_hits = fivep_hits+CNAG_fivep_hits
            threep_hits = threep_hits+CNAG_threep_hits
            new_peaks = new_peaks+CNAG_new_peaks
                                                          
            #Write splice site detection rate and new peak discovery rate to file
            with open("{0}_{1}_detection_rate.txt".format(bedgraph_file.split("/")[-1].split(".")[0], chromosome),"a") as f_handle:
                line = CNAG+"\t"+str(len(CNAG_fivep_hits))+"\t"+str(len(splice_site_dict[CNAG][0]))+"\t"+str(len(CNAG_threep_hits))+"\t"+str(len(splice_site_dict[CNAG][1]))+"\t"+str(len(CNAG_new_peaks))+"\n"
                f_handle.write(line)
            
            with open("{0}_{1}_indexes.txt".format(bedgraph_file.split("/")[-1].split(".")[0], chromosome),"a") as f_handle:
                a = 0
                while a < len(CNAG_fivep_hits): 
                    line_list = [CNAG,transcript_dict[CNAG][3],"known 5' site",str(CNAG_fivep_hits[a]),str(CNAG_fivep_heights[a]),"\n"]
                    line = "\t".join(line_list)
                    f_handle.write(line)
                    a += 1
                b = 0
                while b < len(CNAG_threep_hits):
                    line_list = [CNAG,transcript_dict[CNAG][3],"known 3' site",str(CNAG_threep_hits[b]),str(CNAG_threep_heights[b]),"\n"]
                    line = "\t".join(line_list)
                    f_handle.write(line)
                    b += 1
                c = 0
                while c < len(CNAG_new_peaks):
                    line_list = [CNAG,transcript_dict[CNAG][3],"unknown site",str(CNAG_new_peaks[c]),str(CNAG_new_peaks_heights[c]),"\n"]
                    line = "\t".join(line_list)
                    f_handle.write(line)
                    c += 1
            
    if chromosome != "All":
        print "Chromosome"+chromosome[-1:]
    print "Totals"
    print str(len(fivep_hits))+" out of "+str(fivep_ss_total)+" 5'-splice sites found"
    print str(len(threep_hits))+" out of "+str(threep_ss_total)+" 3'-splice sites found"
    print str(len(new_peaks))+" new peaks"
    

###########################################################################################################
## Tools for recovering sequences in the region around peaks. Peak picking is not sure accurate and has  ##
## off by 1-2 errors. This will take some optimization to land on the actual splice site each time.      ##
###########################################################################################################

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
        return complement(s[::-1])
    
def find_new_peak_sequence(fasta_file, gff3_file, index_file):
    with open("{0}_sequences.txt".format(index_file.split(".")[0]), "w") as fout:
        line = "transcript\t chromosome\t peak type\t coordinate\t peak height\t sequence\t strand\n"
        fout.write(line)
    
    #Read fasta file for chromosome into list
    fasta_dict = {}    
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                chr_num = line[1:6].strip()
            else:
                fasta_dict[chr_num]=line.strip()
    
    #Build transcript dictionary: keys are CNAG, values are [start, end, strand, chromosome]
    transcript_dict = build_transcript_dict(gff3_file)
    
    #Find positions of new peaks
    index_dict = {}
    with open(index_file, "r") as indexes:
        for line in indexes:
            columns = re.split(r'\t', line)
            CNAG = columns[0].strip()
            #if columns[2].strip() == "unknown site":
            if line.startswith("CNAG"):
                chromosome = columns[1]
                position = int(columns[3].strip())
                if transcript_dict[CNAG][2] == "+":
                    sequence = fasta_dict[chromosome][(position-3):(position+7)]
                if transcript_dict[CNAG][2] == "-":
                    sequence = fasta_dict[chromosome][(position-6):(position+4)]
                    sequence = reverse_complement(sequence)
                    
                #Build index dictionary: [transcript, chromosome, peak position, peak height, sequence (-4 to +4)
                index_dict[CNAG] = [CNAG, chromosome, position, columns[4], sequence]
                 
                if sequence[4:6] == "GT":
                    site_class = "5'"
                elif sequence[2:4] == "AG":
                    site_class = "3'"
                else:
                    site_class = "unknown"
                    
                with open("{0}_sequences.txt".format(index_file.split(".")[0]), "a") as fout:
                    line_list = [CNAG, chromosome, columns[2].strip(), str(position), str(columns[4]), sequence, site_class, "\n"]
                    line = "\t".join(line_list)
                    fout.write(line)
                
                 
                    
                    

                
               
                
            
