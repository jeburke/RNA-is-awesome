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
import operator

##################################################################################
## Determines splice site locations from gff3 file. Needs to have "chr" format  ##
##################################################################################

def list_splice_sites(gff3_file, chromosome="All", gene_list=None):
    fin = open(gff3_file,"r")
    transcript_dict = build_transcript_dict(gff3_file)
    #Dictionary will have transcript as key and then a list of 5' splice sites and a list of 3' splice sites as values
    splice_site_dict = {}
    n = 1
    for line in fin:
        columns = re.split(r'\t+', line.strip())
        if len(columns) > 1:
            if columns[2] == "mRNA":
                CNAG = columns[8].strip()
                CNAG = CNAG.split("=")[1]
                CNAG = CNAG.split(";")[0]
                if CNAG.endswith("mRNA"):
                    CNAG = CNAG.split("_")[0]
                if CNAG not in splice_site_dict:
                    splice_site_dict[CNAG] = [[],[],columns[1]]
            
            elif columns[2] == "exon":
                intron_flag=False
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
            
            elif columns[2] == "intron":
                intron_flag=True
                CNAG = columns[8].strip()
                CNAG = CNAG.split("=")[1]
                CNAG = CNAG.split(";")[0]
                if CNAG.endswith("mRNA"):
                    CNAG = CNAG.split("_")[0]
                if CNAG not in splice_site_dict:
                    splice_site_dict[CNAG] = [[],[],columns[1]]
                if gene_list is None:
                    if columns[6] == "+":
                        splice_site_dict[CNAG][0].append(int(columns[3])-2)
                        splice_site_dict[CNAG][1].append(int(columns[4])-1)
                    elif columns[6] == "-":
                        splice_site_dict[CNAG][0].append(int(columns[4]))
                        splice_site_dict[CNAG][1].append(int(columns[3])-1)
                else:
                    if CNAG in gene_list:
                        if columns[6] == "+":
                            splice_site_dict[CNAG][0].append(int(columns[3])-2)
                            splice_site_dict[CNAG][1].append(int(columns[4])-1)
                        elif columns[6] == "-":
                            splice_site_dict[CNAG][0].append(int(columns[4]))
                            splice_site_dict[CNAG][1].append(int(columns[3])-1)
    
    #Trim to entries in gene list
    if gene_list is not None:
        splice_site_dict = dict([(CNAG, splice_site_dict[CNAG]) for CNAG in gene_list if CNAG in splice_site_dict])
    
    #Trim to just one chromosome
    if chromosome != "All":
        splice_site_dict = dict([(CNAG, coords) for CNAG, coords in splice_site_dict.iteritems() if coords[2] == chromosome])
    
    print len(splice_site_dict)
    
    if intron_flag is False:
        for transcript, sites in splice_site_dict.iteritems():
            sites5 = sorted(sites[0])
            sites3 = sorted(sites[1])
            if transcript_dict[transcript][2] == "+":
                sites5 = sites5[:-1]
                sites3 = sites3[1:]
            elif transcript_dict[transcript][2] == "-":
                sites5 = sites5[1:]
                sites3 = sites3[:-1]
            splice_site_dict[transcript] = [sites5, sites3]
    
    if intron_flag is True: 
    #    splice_site_dict = {k: v for k, v in splice_site_dict.items() if len(v[0])>0}
    #    print splice_site_dict
        splice_site_dict = {k: v for k, v in splice_site_dict.items() if k.startswith("Y")}
        
    return (splice_site_dict, intron_flag)

def splice_site_seq(fasta_file, gff3_file):
    splice_site_dict, intron_flag = list_splice_sites(gff3_file)
    #Read fasta file for chromosome into list
    fasta_dict = {}
    letters = ['A','C','G','T','N']
    n = 0
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                n += 1
                chr_num = line[1:-1].strip()
                fasta_dict[chr_num] = str()
            elif line[0] in letters and "chr"+str(n) in fasta_dict:
                fasta_dict["chr"+str(n)] = fasta_dict["chr"+str(n)]+line.strip()

    fasta_dict = collections.OrderedDict(sorted(fasta_dict.items()))
    
    #Build transcript dictionary: keys are CNAG, values are [start, end, strand, chromosome]
    transcript_dict = build_transcript_dict(gff3_file)
    for transcript, sites in splice_site_dict.iteritems():
        for i in range(len(sites[0])):
            #Find sequence surrounding peak
            position_5p = sites[0][i]
            position_3p = sites[1][i]
            chromosome = transcript_dict[transcript][3]
            strand = transcript_dict[transcript][2]
            if strand == "+":
                sequence_5p = fasta_dict[chromosome][(position_5p-10):(position_5p+10)]
                sequence_3p = fasta_dict[chromosome][(position_3p-19):(position_3p+5)]
            elif strand == "-":
                sequence_5p = fasta_dict[chromosome][(position_5p-9):(position_5p+11)]
                sequence_5p = reverse_complement(sequence_5p)
                sequence_3p = fasta_dict[chromosome][(position_3p-4):(position_3p+20)]
                sequence_3p = reverse_complement(sequence_3p)
     
            with open("{0}_5p_ss_seq.txt".format(fasta_file.split("/")[-1].split(".")[0]), "a") as fout:
                #fout.write(strand+"\n")
                fout.write(sequence_5p+"\n")
                
            with open("{0}_3p_ss_seq.txt".format(fasta_file.split("/")[-1].split(".")[0]), "a") as fout:
                #fout.write(strand+"\n")
                fout.write(sequence_3p+"\n")
    
    
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
                    CNAG = CNAG.split("=")[1]
                    CNAG = CNAG.split(";")[0]
                    if CNAG.endswith("mRNA"): CNAG = CNAG.split("_")[0]
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
            
            if bed_chr in transcript_by_chr:
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

def read_CNAGsort_bedgraph(CNAGsorted_bedgraph):
    bedgraph_dict = {}
    with open(CNAGsorted_bedgraph, "r") as CNAGbedgraph:
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
    return bedgraph_dict

def build_piranha_dict(transcript_dict, piranha_output, sort_bedgraph=True, sorted_bedgraph=None, bedgraph_file=None, max_p_value=1.0e-04):
    piranha_dict = {}
    transcript_by_chr = {}
    if sort_bedgraph == True:
        bedgraph_dict = build_bedgraph_dict(transcript_dict, bedgraph_file)
    else:
        bedgraph_dict = read_CNAGsort_bedgraph(sorted_bedgraph)
    transcript_count = 0
    for CNAG, coords in transcript_dict.iteritems():
        if CNAG.endswith("T0"): transcript_count += 1
        elif CNAG.endswith("mRNA"): CNAG = CNAG.split("_")[0]
        chromosome = coords[3]
        piranha_dict[CNAG] = [[],[],[],[]]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(CNAG)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(CNAG)
    
    with open(piranha_output, "r") as fin:
        for line in fin:
            if "nan" not in line:
                columns = re.split(r'\t', line)
                bed_chr = columns[0].strip()
                if columns[5] == "+":
                    bed_position = int(columns[1])
                    strand = "-"
                elif columns[5] == "-":
                    bed_position = int(columns[2])
                    strand = "+"
                bed_peak = float(columns[4])
                p_value = float(columns[6])
            
                if bed_chr in transcript_by_chr:
                    CNAG_list = transcript_by_chr[bed_chr]
                for CNAG in CNAG_list:
                #Dictionary for bedgraph. Values will be [list of genomic positions][# reads starting at that position]
                    if bed_chr == transcript_dict[CNAG][3].strip() and bed_position > transcript_dict[CNAG][0] and bed_position < transcript_dict[CNAG][1] and p_value <= max_p_value and strand == transcript_dict[CNAG][2]:
                        piranha_dict[CNAG][0].append(bed_position)
                        piranha_dict[CNAG][1].append(bed_peak)
                        piranha_dict[CNAG][2].append(p_value)
                        piranha_dict[CNAG][3].append(strand)

    piranha_dict = {k: v for k, v in piranha_dict.items() if len(v[0])>0}
    final_dict = peak_filter(piranha_dict, bedgraph_dict)
    
    with open("{0}_byCNAG.out".format(piranha_output.split(".")[0].split("/")[-1]), "a") as fout:
        fout.write("CNAG\n Position\n Peak height\n P-value\n")
        for CNAG, values in final_dict.iteritems():
            fout.write(CNAG+"\n")
            coord_list = map(str, piranha_dict[CNAG][0])
            coord_line = "\t".join(coord_list)
            fout.write(coord_line+"\n")
            count_list = map(str, piranha_dict[CNAG][1])
            count_line = "\t".join(count_list)
            fout.write(count_line+"\n")
            p_list = map(str, piranha_dict[CNAG][2])
            p_line = "\t".join(p_list)
            fout.write(p_line+"\n")
    
    print transcript_count
    return final_dict

def peak_filter(peak_dict, bedgraph_dict):
    counter1 = 0
    counter2 = 0
    final_dict = {}
    for CNAG, peaks in peak_dict.iteritems():
        final_dict[CNAG] = [[],[],[]]
        for n in range(len(peaks[0])):
            counter1 += 1
            #Find tallest spot in peak in bedgraph
            neighbors = []
            peak_range = range(peaks[0][n]-4,peaks[0][n]+4)
            i=0
            for i in range(len(peak_range)):
                if peak_range[i] in bedgraph_dict[CNAG][0]:
                    bg_index = bedgraph_dict[CNAG][0].index(peak_range[i])
                    neighbors.append((peak_range[i], peaks[1][n], peaks[2][n], bedgraph_dict[CNAG][1][bg_index]))
            if len(neighbors) > 0:
                neighbors = sorted(neighbors, key=operator.itemgetter(3))
                new_peak = neighbors[-1][0]
            
            #pick only peaks with small neighbors and filter by number of reads in peak
                bg_index = bedgraph_dict[CNAG][0].index(new_peak)
                bg_peak = float(bedgraph_dict[CNAG][1][bg_index])
                bg_neigh1 = float(bedgraph_dict[CNAG][1][bg_index-1])
                if bg_index+1 < len(bedgraph_dict[CNAG][1]):
                    bg_neigh2 = float(bedgraph_dict[CNAG][1][bg_index+1])
                if bg_neigh1/bg_peak <= 0.3 and bg_neigh2/bg_peak <= 0.3 and bg_peak > 15 and bg_neigh1/bg_peak >= 0.005:
                    final_dict[CNAG][0].append(new_peak)
                    final_dict[CNAG][1].append(peaks[1][n])
                    final_dict[CNAG][2].append(peaks[2][n])
                    counter2 += 1
                #else: 
                #    print CNAG
                #    print new_peak

    
    print counter1
    print str(counter2)+" peaks"
    return final_dict

###############################################################################################################
## Function for checking the splice sites found by Piranha or the peaks_by_gene function. Takes a dictionary ##
## generated by peaks_by_gene or by build_piranha_dict. Can take a gene list and outputs a tsv file          ##
###############################################################################################################

def check_splice_sites(peak_dict, gff3_file, chromosome="All", gene_list=None, prefix="New"):
    #Create dictionary of known splice sites by transcript
    splice_site_dict, intron_flag = list_splice_sites(gff3_file, chromosome, gene_list=gene_list)
    five_counter = 0
    five_found_counter = 0
    three_counter = 0
    three_found_counter = 0
    new_dict = {}
    new_counter = 0
    counter = 0
    
    for CNAG, splice_sites in splice_site_dict.iteritems():
        if CNAG in peak_dict and (CNAG[-2] != "T" or CNAG.endswith("T0")):
            if CNAG not in new_dict:
                #peak position, peak height, p_value, classification
                new_dict[CNAG]=[[],[],[],[]]
                n = 0
            while n < len(peak_dict[CNAG][0]):
                #new_dict[CNAG][0].append(peak_dict[CNAG][0][n])
                new_dict[CNAG][1].append(peak_dict[CNAG][1][n])
                new_dict[CNAG][2].append(peak_dict[CNAG][2][n])
                peak_range = range(peak_dict[CNAG][0][n]-3, peak_dict[CNAG][0][n]+4)
                five_prime = set(splice_sites[0]).intersection(peak_range)
                three_prime = set(splice_sites[1]).intersection(peak_range)
                if len(five_prime) > 0:
                    five_found_counter += 1
                    five_prime = next(iter(five_prime))
                    new_dict[CNAG][0].append(five_prime)
                    new_dict[CNAG][3].append("5' splice site")
                elif len(three_prime) > 0:
                    three_prime = next(iter(three_prime))
                    new_dict[CNAG][0].append(three_prime)
                    three_found_counter += 1
                    new_dict[CNAG][3].append("3' splice site")
                elif len(five_prime) == 0 and len(three_prime) == 0:
                    new_counter += 1
                    new_dict[CNAG][0].append(peak_dict[CNAG][0][n])
                    new_dict[CNAG][3].append("Unknown")
                n+=1
                

            for five_site in splice_sites[0]:
                site_range = range(five_site-3, five_site+4)
                match = set(site_range).intersection(peak_dict[CNAG][0])
                if len(match) == 0:
                    new_dict[CNAG][0].append(five_site)
                    new_dict[CNAG][1].append(0)
                    new_dict[CNAG][2].append("Not found")
                    new_dict[CNAG][3].append("5' splice site")
                    counter += 1
            for three_site in splice_sites[1]:
                site_range = range(three_site-3, three_site+4)
                match = set(site_range).intersection(peak_dict[CNAG][0])
                if len(match) == 0:
                    new_dict[CNAG][0].append(three_site)
                    new_dict[CNAG][1].append(0)
                    new_dict[CNAG][2].append("Not found")
                    new_dict[CNAG][3].append("3' splice site")
                    counter += 1

            for five_prime in splice_sites[0]:
                five_counter += 1           
            for three_prime in splice_sites[1]:
                three_counter += 1
            
            if intron_flag is False:
                five_counter += -1
                three_counter += -1
            
            
    transcript_counter = 0
    with open("{}_peak_picking.txt".format(prefix), "w") as fout:
        fout.write("Transcript\t Location\t Peak height\t P-value\t Classification\n")
        for CNAG, values in new_dict.iteritems():
            transcript_counter += 1
            n = 0
            while n < len(values[0]):
                line_list = [CNAG, values[0][n], values[1][n], values[2][n], values[3][n], "\n"]
                line_list = map(str, line_list)
                line = "\t".join(line_list)
                fout.write(line)
                n+=1
    
    fivep_percent = float(five_found_counter)/five_counter*100
    threep_percent = float(three_found_counter)/three_counter*100
    print str(five_found_counter)+" out of "+str(five_counter)+" 5' splice sites found, "+str("%0.1f" % fivep_percent)+"%"
    print str(three_found_counter)+" out of "+str(three_counter)+" 3' splice sites found, "+str("%0.1f" % threep_percent)+"%"
    print "...in "+str(transcript_counter)+" transcripts"
    print str(new_counter)+" new peaks found"
    #print counter
    return new_dict

################################################################################################################
## Function to put together peaks from different samples. Grabs all peaks from either sample.                 ##
## Returns dictionary of final peak set                                                                       ##
################################################################################################################

def cat_peaks(dict_list, gff3_file, all_transcripts=False):
    transcript_set = set()
    splice_site_dict, intron_flag = list_splice_sites(gff3_file)
    for peak_dict in dict_list:
        if len(transcript_set) == 0:
            transcript_set = set(peak_dict.keys())
        else:
            if all_transcripts is False:
                new_set = set(peak_dict.keys())
                transcript_set = transcript_set.intersection(new_set)
            elif all_transcripts is True:
                new_set = set(peak_dict.keys())
                transcript_set = transcript_set.union(new_set)
    
    cat_dict = {}
    fivep_total = 0
    threep_total = 0
    fivep_found = 0
    threep_found = 0
    new_found = 0
    
    for transcript in transcript_set:
        fivep_total += len(splice_site_dict[transcript][0])
        threep_total += len(splice_site_dict[transcript][1])
        cat_dict[transcript] = [[],[],[],[]]
        for peak_dict in dict_list:
            if transcript not in peak_dict:
                continue
            #print len(peak_dict)
            i = 0
            for i in range(len(peak_dict[transcript][0])):
                if peak_dict[transcript][0][i] in cat_dict[transcript][0] or peak_dict[transcript][2][i] == "Not found":
                    continue
                else:
                    cat_dict[transcript][0].append(peak_dict[transcript][0][i])
                    cat_dict[transcript][1].append(peak_dict[transcript][1][i])
                    cat_dict[transcript][2].append(peak_dict[transcript][2][i])
                    cat_dict[transcript][3].append(peak_dict[transcript][3][i])
                    if peak_dict[transcript][3][i] == "5' splice site":
                        fivep_found += 1
                    elif peak_dict[transcript][3][i] == "3' splice site":
                        threep_found += 1
                    elif peak_dict[transcript][3][i] == "Unknown":
                        new_found += 1
    
    fivep_percent = float(fivep_found)/fivep_total*100
    threep_percent = float(threep_found)/threep_total*100
    print str(fivep_found)+" out of "+str(fivep_total)+" 5' splice sites found, "+str("%0.1f" % fivep_percent)+"%"
    print str(threep_found)+" out of "+str(threep_total)+" 3' splice sites found, "+str("%0.1f" % threep_percent)+"%"
    print "...in "+str(len(transcript_set))+" transcripts"
    print str(new_found)+" new peaks found"
    #print counter
    return cat_dict
    
def compare_peaks(dict_list, gff3_file, transcript_dict, all_transcripts=False):
    transcript_set = set()
    splice_site_dict, intron_flag = list_splice_sites(gff3_file)
    for peak_dict in dict_list:
        if len(transcript_set) == 0:
            transcript_set = set(peak_dict.keys())
        else:
            if all_transcripts is False:
                new_set = set(peak_dict.keys())
                transcript_set = transcript_set.intersection(new_set)
            elif all_transcripts is True:
                new_set = set(peak_dict.keys())
                transcript_set = transcript_set.union(new_set)
    
    new_dict = {}
    fivep_total = 0
    threep_total = 0
    fivep_found = 0
    threep_found = 0
    new_found = 0

    for transcript in transcript_set:
        peak_set = set()
        fivep_total += len(splice_site_dict[transcript][0])
        threep_total += len(splice_site_dict[transcript][1])
        new_dict[transcript] = [[],[],[],[]]
        chrom = transcript_dict[transcript][3]
        for peak_dict in dict_list:
            peak_list = peak_dict[transcript][0]
            new_set = set((chrom , x) for x in peak_list)
            peak_set = peak_set.union(new_set)

        for peak_dict in dict_list:
            i = 0
            for i in range(len(peak_dict[transcript][0])):
                peak = peak_dict[transcript][0][i]
                height = peak_dict[transcript][1][i]
                p_value = peak_dict[transcript][2][i]
                site_class = peak_dict[transcript][3][i]
                if peak in new_dict[transcript][0] or p_value == "Not found":
                    continue
                elif (chrom, peak) in peak_set:
                    new_dict[transcript][0].append(peak)
                    new_dict[transcript][1].append(height)
                    new_dict[transcript][2].append(p_value)
                    new_dict[transcript][3].append(site_class)
                    if site_class == "5' splice site":
                        fivep_found += 1
                    elif site_class == "3' splice site":
                        threep_found += 1
                    elif site_class == "Unknown":
                        new_found += 1
    
    fivep_percent = float(fivep_found)/fivep_total*100
    threep_percent = float(threep_found)/threep_total*100
    print str(fivep_found)+" out of "+str(fivep_total)+" 5' splice sites found, "+str("%0.1f" % fivep_percent)+"%"
    print str(threep_found)+" out of "+str(threep_total)+" 3' splice sites found, "+str("%0.1f" % threep_percent)+"%"
    print "...in "+str(len(transcript_set))+" transcripts"
    print str(new_found)+" new peaks found"
    #print counter
    return new_dict
    
######################################################################################################
## Pick peaks using PeakUtils package. Can be further refined or substituted with another package   ##
######################################################################################################

def peaks_by_gene(gff3_file, bedgraph_file, chromosome="All", gene_list=None, sort_bedgraph=False, cutoff=1000, show_plots=False):
    #Create dictionary of known splice sites by transcript
    splice_site_dict, intron_flag = list_splice_sites(gff3_file, chromosome, gene_list=gene_list)
    
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
        
        transcript_dict = dict([(CNAG, transcript_dict[CNAG]) for CNAG in gene_list if CNAG in transcript_dict])
    
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
    peak_dict = {}
    CNAG_counter = 0
    final_dict = {}
    counter1 = 0
    counter2 = 0
    for CNAG, values in bedgraph_dict.iteritems():
        bygene_x = numpy.array(values[0])
        bygene_y = numpy.array(values[1])
        
        if numpy.sum(bygene_y) > cutoff:
            peak_dict[CNAG] = [[],[],[]]
            CNAG_counter += 1
            base = peakutils.baseline(bygene_y, 2)
            indexes = peakutils.indexes(bygene_y-base, thres=0.02, min_dist=5)
            
            if gene_list is not None and show_plots is True:
                print CNAG
                pyplot.figure(figsize=(15,6))
                pplot(bygene_x, bygene_y-base, indexes)
            
            #Filter out peaks with less than 15 reads
            a = 0
            bygene_x_filt = []
            bygene_y_filt = []
            while a < len(bygene_x[indexes]):
                if bygene_y[indexes][a] > 15:
                    bygene_x_filt.append(bygene_x[indexes][a])
                    bygene_y_filt.append(bygene_y[indexes][a])
                    if CNAG not in peak_dict: peak_dict[CNAG]=[[],[],[]]
                    peak_dict[CNAG][0].append(bygene_x[indexes][a])
                    peak_dict[CNAG][1].append(bygene_y[indexes][a])
                    peak_dict[CNAG][2].append(0)
                a += 1
                        
            final_dict[CNAG] = [[],[],[]]
            peaks = peak_dict[CNAG]
            for n in range(len(peaks[0])):
                counter1 += 1
                neighbors = []
                peak_range = range(peaks[0][n]-3,peaks[0][n]+3)
                peak_range = [x for x in peak_range if x != peaks[0][n]]
                i=0
                for i in range(len(peaks[0])):
                    if peaks[0][i] in peak_range:
                        neighbors.append((peaks[0][i],peaks[1][i],peaks[2][i]))
                if len(neighbors) > 0:
                    neighbors = sorted(neighbors, key=operator.itemgetter(1))
                    #print neighbors
                    if neighbors[0][0] not in final_dict[CNAG][0]:
                        final_dict[CNAG][0].append(neighbors[0][0])
                        final_dict[CNAG][1].append(neighbors[0][1])
                        final_dict[CNAG][2].append(neighbors[0][2])
                        counter2 += 1
                else:
                    final_dict[CNAG][0].append(peaks[0][n])
                    final_dict[CNAG][1].append(peaks[1][n])
                    final_dict[CNAG][2].append(peaks[2][n])
                    counter2 += 1    
                        
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
                        peak_range = range(bygene_x_filt[i]-3, bygene_x_filt[i]+3)
                        if any(peak in sites[0] for peak in peak_range):
                            CNAG_fivep_hits.append(bygene_x_filt[i])
                            CNAG_fivep_heights.append(bygene_y_filt[i])
                        elif any(peak in sites[1] for peak in peak_range):
                            CNAG_threep_hits.append(bygene_x_filt[i])
                            CNAG_threep_heights.append(bygene_y_filt[i])
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
    print counter1
    print counter2
    return final_dict
    

############################################################################################################
## Tools for recovering sequences in the region around peaks. Peak picking is not super accurate and has  ##
## off by 1-2 errors. Currently looks in range around peak for the splice site                            ##
############################################################################################################

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
        return complement(s[::-1])
    
def find_new_peak_sequence(fasta_file, gff3_file, peak_dict1, peak_dict2=None, prefix="Peak"):
    with open("{0}_sequences.txt".format(prefix), "w") as fout:
        line = "transcript\t chromosome\t peak type\t coordinate\t peak height\t sequence\t Looks like\n"
        fout.write(line)

    #Read fasta file for chromosome into list
    fasta_dict = {}
    letters = ['A','C','G','T','N']
    n = 0
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                n += 1
                chr_num = line[1:-1].strip()
                fasta_dict[chr_num] = str()
            elif line[0] in letters and "chr"+str(n) in fasta_dict:
                fasta_dict["chr"+str(n)] = fasta_dict["chr"+str(n)]+line.strip()

    fasta_dict = collections.OrderedDict(sorted(fasta_dict.items()))
    
    #Build transcript dictionary: keys are CNAG, values are [start, end, strand, chromosome]
    transcript_dict = build_transcript_dict(gff3_file)
    
    #Find positions of new peaks - first file
    index_dict1 = {}
    unknown_counter = 0
    counter5 = 0
    counter3 = 0
    unk5 = 0
    unk3 = 0
    unk = 0
    for CNAG, peaks in peak_dict1.iteritems():
        chromosome = transcript_dict[CNAG][3]
        
        i = 0
        for i in range(len(peaks[0])):
            #Find sequence surrounding peak
            position = peaks[0][i]
            classification = peaks[3][i]
            if transcript_dict[CNAG][2] == "+":
                sequence = fasta_dict[chromosome][(position-2):(position+6)]
                #sequence = fasta_dict[chromosome][(position-10):(position+10)]
            elif transcript_dict[CNAG][2] == "-":
                sequence = fasta_dict[chromosome][(position-5):(position+3)]
                #sequence = fasta_dict[chromosome][(position-10):(position+10)]
                sequence = reverse_complement(sequence)

            #Classify sequences as splice sites
            #if sequence[4:6] == "GT":
            if sequence[3:5] == "GT" or sequence[3:5] == "GC":
                site_class = "5'"
                GUpos = 3
                counter5 += 1
            elif sequence[1:3] == "AG":
                site_class = "3'"
                AGpos = 1
                counter3 += 1
            elif sequence[3:5] != "GT" and sequence[3:5] != "GC" and "GT" in sequence:
                site_class = "5' offset"
                GUpos = sequence.index('GT')
                counter5 += 1
            #elif sequence[2:4] == "AG":
            elif sequence[3:5] != "GT" and sequence[3:5] != "GC" and "AG" in sequence:
                site_class = "3' offset"
                AGpos = sequence.index('AG')
                counter3 += 1
            else:
                site_class = "Unknown"
                unknown_counter += 1
            
            if site_class == "5' offset":
                if transcript_dict[CNAG][2] == "+":
                    sequence = fasta_dict[chromosome][(position+GUpos-5):(position+GUpos+3)]
                elif transcript_dict[CNAG][2] == "-":
                    sequence = fasta_dict[chromosome][(position-GUpos-2):(position-GUpos+6)]
                    sequence = reverse_complement(sequence)
            elif site_class == "3' offset":
                if transcript_dict[CNAG][2] == "+":
                    sequence = fasta_dict[chromosome][(position+AGpos-3):(position+AGpos+5)]
                elif transcript_dict[CNAG][2] == "-":
                    sequence = fasta_dict[chromosome][(position-AGpos-4):(position-AGpos+4)]
                    sequence = reverse_complement(sequence)
                
            if classification == "Unknown" and site_class == "5'": unk5 += 1
            elif classification == "Unknown" and site_class == "3'": unk3 += 1
            elif classification == "Unknown" and site_class == "Unknown": unk += 1

            #Build index dictionary: [transcript, chromosome, known site?, peak position, peak height, p-value, sequence (-4 to +4), classification
            if CNAG in index_dict1:      
                index_dict1[CNAG].append([CNAG, chromosome, peaks[3][i], position, peaks[1][i], peaks[2][i], sequence, site_class])
            else:
                index_dict1[CNAG] = []
                index_dict1[CNAG].append([CNAG, chromosome, peaks[3][i], position, peaks[1][i], peaks[2][i], sequence, site_class])
    
    #Find positions of new peaks - second file   
    if peak_dict2 is not None:
        index_dict2 = {}
        for CNAG, peaks in peak_dict2.iteritems():
            chromosome = transcript_dict[CNAG][3]
            i = 0
            for i in range(len(peaks[0])):
                #Find sequence surrounding peak
                position = peaks[0][i]
                
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
                    site_class = "Unknown"

                #Build index dictionary: [transcript, chromosome, known site?, peak position, peak height, sequence (-4 to +4), classification
                if CNAG in index_dict2:      
                    index_dict2[CNAG].append([CNAG, chromosome, peaks[3][i], position, peaks[1][i], peaks[2][i], sequence, site_class])
                else:
                    index_dict2[CNAG] = []
                    index_dict2[CNAG].append([CNAG, chromosome, peaks[3][i], position, peaks[1][i], peaks[2][i], sequence, site_class])
    
    #Select only peaks in both replicates
        merged_dict = {}
        for CNAG, peaks in index_dict1.iteritems():
            if CNAG in index_dict2:
                merged_dict[CNAG] = []
                for peak in peaks:
                    for peak2 in index_dict2[CNAG]:
                        if peak[:3] == peak2[:3]:
                            merged_dict[CNAG].append(peak)
        index_dict1 = merged_dict
    
    #Make output file
    print "Output location: {0}_sequences.txt".format(prefix)
    
    print "5' splice sites: "+str(counter5)
    print "3' splice sites: "+str(counter3)
    print "Sequences with no splice sites: "+str(unknown_counter)
    
    print "Of unannotated splice sites:"
    print str(unk5)+" 5' splice sites, "+str(unk3)+" 3' splice sites, "+str(unk)+" other sites"
    

    with open("{0}_sequences.txt".format(prefix), "a") as fout:
        for CNAG, peaks in index_dict1.iteritems():
            for peak in peaks:
                peak = map(str, peak)
                line = "\t".join(peak)
                fout.write(line+"\n")
                
#####################################################################################
## Trim bed file reads to just the 5' ends. Helpful for picking peaks in A samples ##
#####################################################################################

def convert_bed_file(bed_file):
    counter = 0
    with open(bed_file, "r") as fin:
        for line in fin:
            counter += 1
            columns = re.split(r'\t', line)
            new_list = []
            new_list.append(columns[0])
            if columns[5].strip() == "+":
                new_list.append(columns[1])
                new_list.append(str(int(columns[1])+1))
            elif columns[5].strip() == "-":
                new_list.append(str(int(columns[2])-1))
                new_list.append(columns[2])
            else:
                print "unknown strand"
                print counter
            new_list = new_list+columns[3:6]
            with open("{0}_5p.bed".format(bed_file.split("/")[-1].split(".")[0]), "a") as fout:
                new_line = "\t".join(new_list)
                fout.write(new_line)

                ###################################################################################################################
## Splits bed file into individual bed file for each transcript - make sure to run this in a new folder because  ##
## it creates the same number of files as transcripts in the organism.                                           ##
###################################################################################################################

def bed_by_CNAG(bed_file, transcript_dict, gene_list=None):
    bed_dict = {}
    transcript_by_chr = {}
    transcript_count = 0
    if gene_list is not None:
        transcript_dict = {key: transcript_dict[key] for key in gene_list}
    
    for CNAG, coords in transcript_dict.iteritems():
        if CNAG.endswith("T0"): transcript_count += 1
        elif CNAG.endswith("mRNA"): CNAG = CNAG.split("_")[0]
        
        chromosome = coords[3]
        bed_dict[CNAG] = []
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(CNAG)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(CNAG)
    
    with open(bed_file, "r") as fin:
        for line in fin:
            columns = re.split(r'\t', line)
            bed_chr = columns[0].strip()
            bed_position1 = int(columns[1])
            bed_position2 = int(columns[2])
            bed_value = columns[3]
            bed_peak = float(columns[4])
            strand = columns[5].strip()
            if strand == "+":
                read_start = int(columns[1])
            elif strand == "-":
                read_start = int(columns[2])
            else: print line
            
            if bed_chr in transcript_by_chr:
                CNAG_list = transcript_by_chr[bed_chr]
            for CNAG in CNAG_list:
                #Dictionary for bedgraph. Values will be [list of genomic positions][reads starting at that position]
                if bed_chr == transcript_dict[CNAG][3].strip() and read_start > transcript_dict[CNAG][0] and read_start < transcript_dict[CNAG][1]:
                    bed_dict[CNAG].append([bed_chr,bed_position1,bed_position2,bed_value,bed_peak,strand,"\n"])
       
    for CNAG, values in bed_dict.iteritems():
        with open("{0}_{1}.bed".format(CNAG, bed_file.split("/")[-1].split(".")[0]), "a") as fout:
            for value in values:
                value = map(str, value)
                line = "\t".join(value)
                fout.write(line)

            
