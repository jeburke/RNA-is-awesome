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
from scipy.optimize import curve_fit
from scipy.misc import factorial
from datetime import datetime

##################################################################################
## Determines splice site locations from gff3 file. Needs to have "chr" format  ##
##################################################################################

def list_splice_sites(gff3_file, chromosome="All", gene_list=None):
    fin = open(gff3_file,"r")
    #Dictionary will have transcript as key and then a list of 5' splice sites and a list of 3' splice sites as values
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

############################################################
## Read bedgraph and sort by transcript into dictionary   ##
############################################################

def build_bedgraph_dict(transcript_dict, bedgraph_file):
    print datetime.now()
    bedgraph_dict = {}
    transcript_by_chr = {}
    for CNAG, coords in transcript_dict.iteritems():
        chromosome = coords[3]
        bedgraph_dict[CNAG] = [[],[]]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(CNAG)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(CNAG)
    
    with open(bedgraph_file, "r") as bedgraph:
        for line in bedgraph:
            columns = re.split(r'\t', line)
            bed_chr = columns[0].strip()
            bed_position = int(columns[1])
            bed_peak = int(columns[3])
            
            CNAG_list = transcript_by_chr[bed_chr]
            for CNAG in CNAG_list:  
                #Dictionary for bedgraph. Values will be [list of genomic positions][reads starting at that position]
                if bed_chr == transcript_dict[CNAG][3].strip() and bed_position > transcript_dict[CNAG][0] and bed_position < transcript_dict[CNAG][1]:
                    bedgraph_dict[CNAG][0].append(bed_position)
                    bedgraph_dict[CNAG][1].append(bed_peak)
   
    with open("{0}_CNAGsort.bedgraph".format(bedgraph_file.split(".")[0]), "a") as fout:
        for CNAG, values in bedgraph_dict.iteritems():
            fout.write(CNAG+"\n")
            coord_list = map(str, bedgraph_dict[CNAG][0])
            coord_line = "\t".join(coord_list)
            fout.write(coord_line+"\n")
            count_list = map(str, bedgraph_dict[CNAG][1])
            count_line = "\t".join(count_list)
            fout.write(count_line+"\n")
            
                    
    bedgraph_dict = collections.OrderedDict(sorted(bedgraph_dict.items()))

    print datetime.now()
    return bedgraph_dict


######################################################################################################
## Pick peaks using PeakUtils package. Can be further refined or substituted with another package   ##
######################################################################################################

def peaks_by_gene(gff3_file, bedgraph_file, chromosome="All", gene_list=None, sort_bedgraph=False, cutoff=1000):
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
    if sort_bedgraph is True:
        print "Please be patient, sorting bedgraph by CNAG"
        bedgraph_dict = build_bedgraph_dict(transcript_dict, bedgraph_file)
    
    #Open file containing sorted bedgraph and build dictionary
    elif sort_bedgraph is False:
        bedgraph_dict = {}
        with open("{0}_CNAGsort.bedgraph".format(bedgraph_file.split(".")[0]), "r") as CNAGbedgraph:
            for k, line in enumerate(CNAGbedgraph):
                line = line.strip()
                if k % 3 == 0:
                    CNAG = line
                    bedgraph_dict[CNAG] = []
                if k % 3 == 1:
                    coord_list = line.split("\t")
                    while '' in coord_list:
                        coord_list.remove('')
                    coord_list = map(int, coord_list)
                    bedgraph_dict[CNAG].append(coord_list)
                if k % 3 == 2:
                    count_list = line.split("\t")
                    while '' in count_list:
                        count_list.remove('')
                    count_list = map(int, count_list)
                    bedgraph_dict[CNAG].append(count_list)
        
    
    #Create numpy arrays to pick peaks from bedgraph dictionary
    CNAG_counter = 0
    for CNAG, values in bedgraph_dict.iteritems():
        bygene_x = numpy.array(values[0])
        bygene_y = numpy.array(values[1])
        
        if numpy.sum(bygene_y) > cutoff:
            CNAG_counter += 1
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
    print str(CNAG_counter)+" transcripts above cutoff"
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
    
def find_new_peak_sequence(fasta_file, gff3_file, index_file1, index_file2):
    with open("{0}_sequences.txt".format(index_file1.split(".")[0]), "w") as fout:
        line = "transcript\t chromosome\t peak type\t coordinate\t peak height\t sequence\t Looks like\n"
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
    
    #Find positions of new peaks - first file
    index_dict1 = {}
    i = 0
    with open(index_file1, "r") as indexes:
        for line in indexes:
            if i == 0:
                i += 1
            else:
                columns = re.split(r'\t', line)
                CNAG = columns[0].strip()
                chromosome = columns[1]
                position = int(columns[3].strip())
            
                #Find sequence surrounding peak
                if transcript_dict[CNAG][2] == "+":
                    sequence = fasta_dict[chromosome][(position-3):(position+7)]
                elif transcript_dict[CNAG][2] == "-":
                    sequence = fasta_dict[chromosome][(position-6):(position+4)]
                    sequence = reverse_complement(sequence)
                        
                #Classify sequences as splice sites
                if sequence[4:6] == "GT":
                    site_class = "5'"
                elif sequence[2:4] == "AG":
                    site_class = "3'"
                else:
                    site_class = "unknown"

                #Build index dictionary: [transcript, chromosome, known site?, peak position, peak height, sequence (-4 to +4), classification
                if CNAG in index_dict1:      
                    index_dict1[CNAG].append([CNAG, chromosome, columns[2], position, columns[4], sequence, site_class])
                else:
                    index_dict1[CNAG] = []
                    index_dict1[CNAG].append([CNAG, chromosome, columns[2], position, columns[4], sequence, site_class])
    
    #Find positions of new peaks - second file   
    index_dict2 = {}
    j = 0
    with open(index_file2, "r") as indexes:
        for line in indexes:
            if j == 0:
                j += 1
            else:
                columns = re.split(r'\t', line)
                CNAG = columns[0].strip()
                chromosome = columns[1]
                position = int(columns[3].strip())
            
                #Find sequence surrounding peak
                if transcript_dict[CNAG][2] == "+":
                    sequence = fasta_dict[chromosome][(position-3):(position+7)]
                elif transcript_dict[CNAG][2] == "-":
                    sequence = fasta_dict[chromosome][(position-6):(position+4)]
                    sequence = reverse_complement(sequence)
                        
                #Classify sequences as splice sites
                if sequence[4:6] == "GT":
                    site_class = "5'"
                elif sequence[2:4] == "AG":
                    site_class = "3'"
                else:
                    site_class = "unknown"

                #Build index dictionary: [transcript, chromosome, known site?, peak position, peak height, sequence (-4 to +4), classification
                if CNAG in index_dict2:      
                    index_dict2[CNAG].append([CNAG, chromosome, columns[2], position, columns[4], sequence, site_class])
                else:
                    index_dict2[CNAG] = []
                    index_dict2[CNAG].append([CNAG, chromosome, columns[2], position, columns[4], sequence, site_class])
    
    #Select only peaks in both replicates
    merged_dict = {}
    for CNAG, peaks in index_dict1.iteritems():
        if CNAG in index_dict2:
            merged_dict[CNAG] = []
            for peak in peaks:
                for peak2 in index_dict2[CNAG]:
                    if peak[:3] == peak2[:3]:
                        merged_dict[CNAG].append(peak)
    
    #Make output file
    print "Output location: {0}_sequences.txt".format(index_file1.split(".")[0])
    for CNAG, peaks in merged_dict.iteritems():
        with open("{0}_sequences.txt".format(index_file1.split(".")[0]), "a") as fout:
            for peak in peaks:
                peak = map(str, peak)
                line = "\t".join(peak)
                fout.write(line+"\n")
                
###########################################################################################
## Peak picking option #2 - remove known peaks, fit remaining to Poisson distribution    ##
## then find outliers                                                                    ##
###########################################################################################

def poisson_peaks(gff3_file, bedgraph_file, chromosome="All", gene_list=None):
    #Create dictionary of known splice sites by transcript
    splice_site_dict = list_splice_sites(gff3_file, chromosome, gene_list=gene_list)

    #Build transcript dictionary
    transcript_dict = build_transcript_dict(gff3_file)

    if gene_list is not None:
        transcript_dict = dict([(CNAG, transcript_dict[CNAG]) for CNAG in gene_list])

    #Read bedgraph and sort by transcript into dictionary
    bedgraph_dict = build_bedgraph_dict(transcript_dict, bedgraph_file)

    #Remove known peaks from bedgraph
    bg_bedgraph_dict = {}
    for CNAG, values in bedgraph_dict.iteritems():
        print CNAG
        bg_bedgraph_dict[CNAG] = [[],[]]
        n = 0
        while n < len(values[0]):
            if values[0][n] in splice_site_dict[CNAG][0] or values[0][n] in splice_site_dict[CNAG][1]:
                n += 1
            else:
                bg_bedgraph_dict[CNAG][0].append(values[0][n])
                bg_bedgraph_dict[CNAG][1].append(values[1][n])
                n += 1

    value_list = []
    for CNAG, values in bedgraph_dict.iteritems():
        x = numpy.array(values[0])
        y = numpy.array(values[1])
        value_list = value_list + values[1]

        #if gene_list is not None:
            #pyplot.figure(figsize=(15,6))
            #pyplot.plot(x,y)


    background_list = []
    for CNAG, values in bg_bedgraph_dict.iteritems():
        bg_x = numpy.array(values[0])
        bg_y = numpy.array(values[1])
        background_list = background_list + values[1]

        #if gene_list is not None:
            #pyplot.figure(figsize=(15,6))
            #pyplot.plot(bg_x, bg_y)

    #Remove zeros from background data
    print len(background_list)
    background_list = list(filter(lambda x: x>0, background_list))
    print len(background_list)
    
    #Poisson fit to background data (minus known splice sites)

    entries, bin_edges, patches = pyplot.hist(background_list, bins=50, range=[0,200], normed=True)
    bin_middles = 0.5*(bin_edges[1:] + bin_edges[:-1])

    #Define function that fits data and return values
    def poisson(k, lamb):
        return(lamb**k/factorial(k) * numpy.exp(-lamb))

    parameters, cov_matrix = curve_fit(poisson, bin_middles, entries)
    print parameters

    x_plot = numpy.linspace(0, 500, 501)
    pyplot.plot(x_plot, poisson(x_plot, *parameters), 'r-', lw=2)
    pyplot.xlim([0,50])
    pyplot.show()

    #Define residuals based on fit
    outlier_list = []
    lamb = parameters[0]
    print "Maximum peak height:"
    print max(background_list)
    print "Median peak height:"
    print numpy.median(numpy.array(background_list))
    for x in x_plot:
        if x > max(background_list)/2:
            expected = lamb**x/factorial(x) * numpy.exp(-lamb)
            observed = 0
            for value in background_list:
                if value == x:
                    observed += 1
            if observed > expected:
                outlier_list.append(x)
    print outlier_list

    



