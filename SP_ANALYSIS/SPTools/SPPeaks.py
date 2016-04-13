__author__ = 'jordanburke'

import pandas
import math
import random
import re
import numpy
import peakutils
from peakutils.plot import plot as pplot
from matplotlib import pyplot

##################################################################################
## Determines splice site locations from gff3 file. Needs to have "chr" format  ##
##################################################################################

def list_splice_sites(gff3_file, chromosome, gene_list=None, sites="All"):
    fin = open(gff3_file,"r")
    fivePsites = []
    threePsites = []
    for line in fin:
        if line.startswith("chr7"):
            columns = re.split(r'\t+', line)
            CNAG = columns[8].strip()
            CNAG = CNAG[-12:]      
            if gene_list is None:
                if columns[2] == "exon":
                    if columns[6] == "+":
                        fivePsites.append(int(columns[4])-1)
                        threePsites.append(int(columns[3])-2)
                    if columns[6] == "-":
                        fivePsites.append(int(columns[3])-1)
                        threePsites.append(int(columns[4]))
            else:
                if CNAG in gene_list and columns[2] == "exon":
                    if columns[6] == "+":
                        fivePsites.append(int(columns[4])-1)
                        threePsites.append(int(columns[3])-2)
                    if columns[6] == "-":
                        fivePsites.append(int(columns[3])-1)
                        threePsites.append(int(columns[4]))
                    
    print "Five prime splice sites: " +str(len(fivePsites))
    print "Three prime splice sites: " +str(len(threePsites))
    allsites = fivePsites+threePsites
    if sites == "All":
        return allsites
    elif sites == "3p":
        return threePsites
    elif sites == "5p":
        return fivePsites
    
######################################################################################################
## Pick peaks using PeakUtils package. Can be further refined or substituted with another package   ##
######################################################################################################

def peaks_by_gene(gff3_file, bedgraph_file, chromosome, gene_list=None):
    gff3 = open(gff3_file,"r")
    transcript_dict = {}
    splice_sites = list_splice_sites(gff3_file, chromosome, gene_list=gene_list)
    for line in gff3:
        if line.startswith(chromosome):
            columns = re.split(r'\t+', line)
            if columns[2] == "mRNA":
                CNAG = columns[8]
                CNAG = CNAG[3:15]
                transcript_dict[CNAG] = [int(columns[3]), int(columns[4]), columns[6]]
    
    if gene_list is not None:
        transcript_dict = dict([(CNAG, transcript_dict[CNAG]) for CNAG in gene_list])

    hits = []
    misses = []
    counter = 0
    for CNAG, coords in transcript_dict.iteritems():
        bedgraph = open(bedgraph_file,"r")
        bygene_x = []
        bygene_y = []
        for line in bedgraph:
            if line.startswith(chromosome):
                columns = re.split(r'\t+', line)
                if int(columns[1]) > coords[0] and int(columns[1]) < coords[1]:
                    bygene_x.append(int(columns[1]))
                    bygene_y.append(int(columns[3]))
        bygene_x = numpy.array(bygene_x)
        bygene_y = numpy.array(bygene_y)
        with open("test_arrays.txt","a") as f_handle:
            numpy.savetxt(f_handle, bygene_x, fmt='%i')
        

        if numpy.sum(bygene_y) > 100:
            indexes = peakutils.indexes(bygene_y, thres=0.05, min_dist=1)
            
            if gene_list is not None:
                print CNAG
                print(bygene_x[indexes])
                pyplot.figure(figsize=(15,6))
                pplot(bygene_x, bygene_y, indexes)
            
            
            with open("{0}_indexes.txt".format(chromosome),"a") as f_handle:
                numpy.savetxt(f_handle, bygene_x[indexes], fmt='%i')
        

            for site in splice_sites:
                if site in bygene_x[indexes]:
                    hits.append(site)
                else:
                    misses.append(site)
            
            counter += 1
    
    print counter
    print "Of "+str(len(splice_sites))+" sites:"
    print str(len(hits))+" hits "#+str(len((allhits)/len(allsites))*100)+"%" 
    #print str(len(allmisses))+" misses," +str(len(allmisses)/len(allsites)*100)+"%"