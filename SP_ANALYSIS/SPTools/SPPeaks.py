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
## Pick peaks using PeakUtils package. Can be further refined or substituted with another package   ##
######################################################################################################

def peaks_by_gene(gff3_file, bedgraph_file, chromosome="All", gene_list=None):
    gff3 = open(gff3_file,"r")
    transcript_dict = {}
    
    #Create dictionary of known splice sites by transcript
    splice_site_dict = list_splice_sites(gff3_file, chromosome, gene_list=gene_list)
    
    #Initiate output files
    with open("{0}_indexes.txt".format(chromosome),"w") as f_handle:
        line = "transcript\t"+chromosome+" coordinate\t"+"peak size\n"
        f_handle.write(line)
    with open("{0}_detection_rate.txt".format(chromosome),"w") as f_handle:
        line = "transcript\t 5' splice site hits\t 5' splice sites \t 3' splice site hits \t 3' splice sites \t new peaks \n"
        f_handle.write(line)
   
    bedgraph_dict = {}
    for line in gff3:
        columns = re.split(r'\t+', line)
        if len(columns) > 1:
            if columns[2] == "mRNA":
                CNAG = columns[8]
                CNAG = CNAG[3:15]
                #Transcript dictionary: keys are CNAG, values are [start, end, strand, chromosome]
                transcript_dict[CNAG] = [int(columns[3]), int(columns[4]), columns[6], columns[0]]
                #Empty dictionary for bedgraph. Values will be [list of genomic positions][reads starting at that position]
                bedgraph_dict[CNAG] = [[],[]]
                
    trasncript_dict = collections.OrderedDict(sorted(transcript_dict.items()))
    bedgraph_dict = collections.OrderedDict(sorted(bedgraph_dict.items()))
    
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
    with open(bedgraph_file, "r") as bedgraph:
        for line in bedgraph:
            columns = re.split(r'\t', line)
            for CNAG, coords in transcript_dict.iteritems():
                 if int(columns[1]) > coords[0] and int(columns[1]) < coords[1]:
                    bedgraph_dict[CNAG][0].append(int(columns[1]))
                    bedgraph_dict[CNAG][1].append(int(columns[3]))
    
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
                
            #Write peaks and peak heights to file
            with open("{0}_{1}_indexes.txt".format(bedgraph_file.split("/")[-1].split(".")[0], chromosome),"a") as f_handle:
                i = 0
                while i < len(bygene_x_filt):
                    line = CNAG+"\t"+str(bygene_x_filt[i])+"\t"+str(bygene_y_filt[i])+"\n"
                    f_handle.write(line)
                    i += 1
                                                          
            #Compare peaks to known splice sites from splice site dictionary                                             
            CNAG_fivep_hits = []
            CNAG_threep_hits = []
            CNAG_new_peaks = []
                                                          
            for dict_CNAG, sites in splice_site_dict.iteritems():
                if dict_CNAG == CNAG:
                    fivep_ss_total += len(sites[0])
                    threep_ss_total += len(sites[1])
                    for peak in bygene_x_filt:
                        if peak in sites[0]:
                            CNAG_fivep_hits.append(peak)
                        elif peak in sites[1]:
                            CNAG_threep_hits.append(peak)
                        else:
                            CNAG_new_peaks.append(peak)
                                                          
            #Add number for this transcript to existing totals
            fivep_hits = fivep_hits+CNAG_fivep_hits
            threep_hits = threep_hits+CNAG_threep_hits
            new_peaks = new_peaks+CNAG_new_peaks
                                                          
            #Write splice site detection rate and new peak discovery rate to file
            with open("{0}_{1}_detection_rate.txt".format(bedgraph_file.split("/")[-1].split(".")[0], chromosome),"a") as f_handle:
                line = CNAG+"\t"+str(len(CNAG_fivep_hits))+"\t"+str(len(splice_site_dict[CNAG][0]))+"\t"+str(len(CNAG_threep_hits))+"\t"+str(len(splice_site_dict[CNAG][1]))+"\t"+str(len(CNAG_new_peaks))+"\n"
                f_handle.write(line)
    
    
    if chromosome != "All":
        print "Chromosome"+chromosome[-1:]
    print "Totals"
    print str(len(fivep_hits))+" out of "+str(fivep_ss_total)+" 5'-splice sites found"
    print str(len(threep_hits))+" out of "+str(threep_ss_total)+" 3'-splice sites found"
    print str(len(new_peaks))+" new peaks"
    
#def find_new_peak_sequence(fasta_file, chromosome, index_file):
    #Read fasta file for chromosome into list
    
    #
    