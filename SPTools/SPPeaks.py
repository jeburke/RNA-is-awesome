__author__ = 'jordanburke'

import sys
import pandas as pd
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
from scipy import stats
from datetime import datetime
import operator
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/SPTools')
import SPTables
import SPPlots
import SPScores
import SPJunctions
from math import log
from matplotlib import pyplot as plt


######################################################################################################
## Build a transcript dictionary with all information from a gff3 file                              ##
######################################################################################################

def build_transcript_dict(gff3_file, organism=None):
    rom_lat = {'I':'chr1','II':'chr2','III':'chr3','IV':'chr4','V':'chr5','VI':'chr6','VII':'chr7','VIII':'chr8','IX':'chr9','X':'chr10','XI':'chr11','XII':'chr12','XIII':'chr13','XIV':'chr14','XV':'chr15','XVI':'chr16','MT':'MT'}
    with open(gff3_file,"r") as gff3:
        transcript_dict = {}
        for line in gff3:
            columns = re.split(r'\t+', line)
            
            if organism == 'pombe' and len(columns) > 1:
                chr_rom = columns[0]
                chrom = rom_lat[chr_rom]
                transcript_types = ['transcript','pseudogene','rRNA','snoRNA','tRNA','snRNA']
                if columns[2] in transcript_types:
                    if columns[8].split(':')[0].split('=')[1] == 'gene': continue
                    transcript = columns[8].split(';')[0].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], chrom, [], []]
                elif columns[2] == 'CDS':
                    transcript = columns[8].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
                        
            if len(columns) == 9 and organism is None:
                if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                    transcript = columns[8]
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"): transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    #Transcript dictionary: keys are transcript, values are [start, end, strand, chromosome, CDS start, CDS end]
                    
                    chrom = columns[0]
                    if chrom in rom_lat: chrom = rom_lat[chrom]
                    elif chrom not in rom_lat.keys():
                        if chrom[3:] in rom_lat: chrom = rom_lat[chrom[3:]]
                            
                    if transcript not in transcript_dict:
                        transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], chrom, [], []]
                    else:
                        transcript_dict[transcript][0] = int(columns[3])
                        transcript_dict[transcript][1] = int(columns[4])
                        transcript_dict[transcript][2] = columns[6]
                        transcript_dict[transcript][3] = columns[0]
                elif columns[2] == "CDS":
                    transcript = columns[8].split("=")[1].split(".")[0].split(';')[0]
                    if 'mRNA' in transcript: transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    if transcript not in transcript_dict:
                        strand = columns[6]
                        chrom = columns[0]
                        transcript_dict[transcript] = [0,0,strand,chrom,[],[]]
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
                    
    for tx in transcript_dict:
        if transcript_dict[tx][0] == 0:
            transcript_dict[tx][0] = transcript_dict[tx][4][0]
            transcript_dict[tx][1] = transcript_dict[tx][5][0]
    
    transcript_dict = collections.OrderedDict(sorted(transcript_dict.items()))
    return transcript_dict


##################################################################################
## Determines splice site locations from gff3 file. Needs to have "chr" format  ##
##################################################################################

def list_splice_sites(gff3_file, chromosome="All", gene_list=None, organism=None):
    fin = open(gff3_file,"r")
    transcript_dict = build_transcript_dict(gff3_file, organism=organism)
    #Dictionary will have transcript as key and then a list of 5' splice sites and a list of 3' splice sites as values
    splice_site_dict = {}
    n = 1
    
    #Read gff3 file and find all exon entries
    for line in fin:
        columns = re.split(r'\t+', line.strip())
        
        if organism == 'pombe' and len(columns)>1:
            intron_flag = False
            if columns[2] == 'exon':
                chr_rom = columns[0]
                rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
                chrom = rom_lat[chr_rom]
                transcript = columns[8].split(':')[0].split('=')[1]
                
                if transcript not in splice_site_dict:
                    splice_site_dict[transcript] = [[],[],chrom]
                if columns[6] == "+":
                    splice_site_dict[transcript][0].append(int(columns[4])-1)
                    splice_site_dict[transcript][1].append(int(columns[3])-2)
                elif columns[6] == "-":
                    splice_site_dict[transcript][0].append(int(columns[3])-1)
                    splice_site_dict[transcript][1].append(int(columns[4]))
        
        if len(columns) > 1 and organism is None:
            if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                transcript = columns[8].strip()
                transcript = transcript.split("=")[1]
                transcript = transcript.split(";")[0]
                if transcript.endswith("mRNA"):
                    transcript = transcript.split("_")[0]
                if transcript[-2] != 'T':
                    transcript = transcript+'T0'
                if transcript not in splice_site_dict:
                    splice_site_dict[transcript] = [[],[],columns[0]]
            
            elif columns[2] == "exon":
                intron_flag=False
                transcript = columns[8].strip()
                transcript = transcript[-12:]
                if transcript[-2] != 'T':
                    transcript = transcript+'T0'
                if gene_list is None:
                    if columns[6] == "+":
                        splice_site_dict[transcript][0].append(int(columns[4])-1)
                        splice_site_dict[transcript][1].append(int(columns[3])-2)
                    if columns[6] == "-":
                        splice_site_dict[transcript][0].append(int(columns[3])-1)
                        splice_site_dict[transcript][1].append(int(columns[4]))
                else:
                    if transcript in gene_list:
                        if columns[6] == "+":
                            splice_site_dict[transcript][0].append(int(columns[4])-1)
                            splice_site_dict[transcript][1].append(int(columns[3])-2)
                        if columns[6] == "-":
                            splice_site_dict[transcript][0].append(int(columns[3])-1)
                            splice_site_dict[transcript][1].append(int(columns[4]))
            
            #For organisms where introns are annotated instead of exons (e.g. S. cerevisiae)
            elif "intron" in columns[2]: 
                intron_flag=True
                transcript = columns[8].strip()
                transcript = transcript.split("=")[1]
                transcript = transcript.split(";")[0]
                if transcript.endswith("mRNA"):
                    transcript = transcript.split("_")[0]
                if transcript[-2] != 'T':
                    transcript = transcript+'T0'
                if transcript not in splice_site_dict:
                    splice_site_dict[transcript] = [[],[],columns[0]]
                if gene_list is None:
                    if columns[6] == "+":
                        splice_site_dict[transcript][0].append(int(columns[3])-1)
                        splice_site_dict[transcript][1].append(int(columns[4]))
                    elif columns[6] == "-":
                        splice_site_dict[transcript][0].append(int(columns[4]))
                        splice_site_dict[transcript][1].append(int(columns[3])-1)
                else:
                    if transcript in gene_list:
                        if columns[6] == "+":
                            splice_site_dict[transcript][0].append(int(columns[3])-2)
                            splice_site_dict[transcript][1].append(int(columns[4])-1)
                        elif columns[6] == "-":
                            splice_site_dict[transcript][0].append(int(columns[4]))
                            splice_site_dict[transcript][1].append(int(columns[3])-1)
    
    
    #Trim to entries in gene list
    if gene_list is not None:
        splice_site_dict = {transcript: splice_site_dict[transcript] for transcript in gene_list}
        #new_dict = {}
        #for transcript, sites in splice_site_dict.iteritems():
        #    if transcript in gene_list:
        #        print transcript
        #        new_dict[transcript] = sites
        #splice_site_dict = new_dict
    #Trim to just one chromosome
    if chromosome != "All":
        splice_site_dict = dict([(transcript, coords) for transcript, coords in splice_site_dict.iteritems() if coords[2] == chromosome])
    
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
    
    #if intron_flag is True: 
    #    splice_site_dict = {k: v for k, v in splice_site_dict.items() if k.startswith("Y")}
        
    return (splice_site_dict, intron_flag)

def splice_site_seq(fasta_file, gff3_file, gene_list=None):
    splice_site_dict, intron_flag = list_splice_sites(gff3_file, gene_list=gene_list)
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
    
    #Build transcript dictionary: keys are transcript, values are [start, end, strand, chromosome]
    transcript_dict = build_transcript_dict(gff3_file)
    for transcript, sites in splice_site_dict.iteritems():
        for i in range(len(sites[0])):
            #Find sequence surrounding splice site
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
    

############################################################
## Read bedgraph and sort by transcript into dictionary   ##
############################################################

def build_bedgraph_dict(transcript_dict, bedgraph_file):
    print datetime.now()
    bedgraph_dict = {}
    transcript_by_chr = {}
    for transcript, coords in transcript_dict.iteritems():
        chromosome = coords[3]
        bedgraph_dict[transcript] = [[],[]]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(transcript)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(transcript)
    
    with open(bedgraph_file, "r") as bedgraph:
        for line in bedgraph:
            columns = re.split(r'\t', line)
            bed_chr = columns[0].strip()
            rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
            if bed_chr in rom_lat:
                bed_chr = rom_lat[bed_chr]
            bed_position = int(columns[1])
            bed_peak = float(columns[3])
            
            if bed_chr in transcript_by_chr:
                transcript_list = transcript_by_chr[bed_chr]
            for transcript in transcript_list:  
                #Dictionary for bedgraph. Values will be [list of genomic positions][reads starting at that position]
                if bed_chr == transcript_dict[transcript][3].strip() and bed_position > transcript_dict[transcript][0] and bed_position < transcript_dict[transcript][1]:
                    bedgraph_dict[transcript][0].append(bed_position)
                    bedgraph_dict[transcript][1].append(bed_peak)
   
    with open("{0}_CNAGsort.bedgraph".format(bedgraph_file.split("/")[-1].split(".")[0]), "a") as fout:
        for transcript, values in bedgraph_dict.iteritems():
            fout.write(transcript+"\n")
            coord_list = map(str, bedgraph_dict[transcript][0])
            coord_line = "\t".join(coord_list)
            fout.write(coord_line+"\n")
            count_list = map(str, bedgraph_dict[transcript][1])
            count_line = "\t".join(count_list)
            fout.write(count_line+"\n")
                                
    bedgraph_dict = collections.OrderedDict(sorted(bedgraph_dict.items()))

    print datetime.now()
    #print bedgraph_dict
    #return bedgraph_dict

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

def read_CNAGsort_bedgraph2(bedgraph_dict_output, transcript_dict, organism=None):
    bg_dict = {}
    count = 0
    dtype = [('coord', int), ('height', float)]
    with open(bedgraph_dict_output,'r') as f:
        n = 0
        for line in f:
            n += 1
            if len(line) > 1:
                
                #Read the transcript line
                if n%3 == 1:
                    tx = line.strip()
                    count += 1
                    #print count
                    if tx[-2] != 'T' and organism != 'pombe':
                        tx = tx+'T0'
                
                #Read the coordinate line
                elif n%3 == 2:
                    coords  = map(int, line.strip().split('\t'))
                
                #Read the values line
                elif n%3 == 0:
                    heights = map(float, line.strip().split('\t'))

                    if tx not in transcript_dict:
                        pass
                    else:
                        all_coords = set(range(min(coords),max(coords)))
                        missing = all_coords.difference(coords)
                        coords = coords + list(missing)

                        #Fill in missing coordinates with zeros
                        zero_fill = [0]*len(missing)
                        heights = heights + zero_fill

                        #Create a pandas series with all coordinates and sort so zeros are inserted appropriately
                        entry = pd.Series(heights, index=coords)
                        entry.sort_index(inplace=True)

                        selected_range = range(transcript_dict[tx][0],transcript_dict[tx][1])
                        entry = entry[entry.index.isin(selected_range)]
                        
                        bg_dict[tx] = entry
    return bg_dict

def build_piranha_dict(transcript_dict, piranha_output, sort_bedgraph=True, sorted_bedgraph=None, bedgraph_file=None, max_p_value=1.0e-04):
    piranha_dict = {}
    transcript_by_chr = {}
    
    #Make dictionary with snRNAs and tRNAs
    nc_dict = {}
    for transcript, values in transcript_dict.iteritems():
        if transcript.startswith("sn") or transcript.startswith("t") or transcript == "NME1":
            nc_dict[transcript]=values
            
    #Sort bedgraph file by CNAG or (to save time) read in already sorted file
    if sort_bedgraph == True:
        bedgraph_dict = build_bedgraph_dict(transcript_dict, bedgraph_file)
    else:
        bedgraph_dict = read_CNAGsort_bedgraph(sorted_bedgraph)
    transcript_count = 0
    
    #Build dictionary with transcripts in each chromosome - this speeds up processing by limiting where the program needs to look for each peak
    for CNAG, coords in transcript_dict.iteritems():
        if CNAG.endswith("T0") or CNAG.endswith('.1'): transcript_count += 1
        elif CNAG.endswith("mRNA"): CNAG = CNAG.split("_")[0]
        chromosome = coords[3]
        rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
        if chromosome in rom_lat:
            chromosome = rom_lat[chromosome]
        piranha_dict[CNAG] = [[],[],[],[]]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(CNAG)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(CNAG)
    
    #Read piranha output into dicionary. Values are [genomic positions][# reads at each position]
    #Also filters by p_value and makes sure that the strand is correct based on the transcript
    with open(piranha_output, "r") as fin:
        for line in fin:
            nc_flag = False
            if line.startswith('A'): continue
            if "nan" not in line:
                columns = re.split(r'\t', line)
                bed_chr = columns[0].strip()
                if bed_chr in rom_lat:
                    bed_chr = rom_lat[bed_chr]
                if columns[5] == "+":
                    bed_position = int(columns[1])
                    strand = "-"
                elif columns[5] == "-":
                    bed_position = int(columns[2])
                    strand = "+"
                bed_peak = float(columns[4])
                p_value = float(columns[6])
            
                #Get list of transcripts for current chromosome
                if bed_chr in transcript_by_chr:
                    CNAG_list = transcript_by_chr[bed_chr]
                
                #Check if peak is in a snoRNA or tRNA
                for nc_RNA, coords in nc_dict.iteritems():
                    if bed_position < coords[1]+5 and bed_position > coords[0]-5:
                        nc_flag = True   
                
                for CNAG in CNAG_list:
                    if bed_chr == transcript_dict[CNAG][3].strip() and bed_position > transcript_dict[CNAG][0] and bed_position < transcript_dict[CNAG][1] and p_value <= max_p_value and strand == transcript_dict[CNAG][2] and nc_flag == False and bed_chr != 'MT':     
                        piranha_dict[CNAG][0].append(bed_position)
                        piranha_dict[CNAG][1].append(bed_peak)
                        piranha_dict[CNAG][2].append(p_value)
                        piranha_dict[CNAG][3].append(strand)

    #Filter piranha dict to get rid of transcripts with no peaks
    piranha_dict = {k: v for k, v in piranha_dict.items() if len(v[0])>0}
    
    #Apply other filters (see function) including minimum number of reads and height compared to neighbors
    final_dict = peak_filter(piranha_dict, bedgraph_dict)
    

    with open("{0}_byCNAG.out".format(piranha_output.split("/")[-1].split(".")[0]), "a") as fout:
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
                if bg_neigh1/bg_peak <= 0.3 and bg_neigh2/bg_peak <= 0.3 and bg_peak >= 10 and bg_neigh1/bg_peak >= 0.001:
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

def check_splice_sites(peak_dict, gff3_file, chromosome="All", gene_list=None, prefix="New", organism=None):
    #Create dictionary of known splice sites by transcript
    splice_site_dict, intron_flag = list_splice_sites(gff3_file, chromosome, gene_list=gene_list, organism=organism)
    five_counter = 0
    five_found_counter = 0
    three_counter = 0
    three_found_counter = 0
    new_peak_dict = {}
    new_dict = {}
    new_counter = 0
    counter = 0
    
    #Collapse isoforms for each gene
    for transcript, splice_sites in splice_site_dict.iteritems():
        if transcript in peak_dict:
            gene = transcript[:-2]
            if gene not in new_peak_dict:
                peak_list = zip(peak_dict[transcript][0], peak_dict[transcript][1])
                peak_set = set(peak_list)
                five_set = set(splice_sites[0])
                three_set = set(splice_sites[1])
            else:
                peak_list = zip(peak_dict[transcript][0], peak_dict[transcript][1])
                peak_set = new_peak_dict[gene][0].union(peak_list)
                five_set = new_peak_dict[gene][1].union(splice_sites[0])
                three_set = new_peak_dict[gene][2].union(splice_sites[1])
            new_peak_dict[gene]=[peak_set,five_set,three_set]
            
    #Build new dictionary that classifies peaks and lists undetected splice sites
    for gene, sets in new_peak_dict.iteritems():
        #peak position, peak height, classification
        new_dict[gene] = [[],[],[]]
        n=0
        peak_tuples = list(new_peak_dict[gene][0])
        for n in range(len(peak_tuples)):
            #new_dict[gene][0].append(new_peak_dict[gene][0][n][0])
            new_dict[gene][1].append(peak_tuples[n][1])
                
            #Create a range surrounding the peak and check if there is an annotated splice site in that range
            peak_range = range(peak_tuples[n][0]-3, peak_tuples[n][0]+4)
            five_prime = new_peak_dict[gene][1].intersection(peak_range)
            three_prime = new_peak_dict[gene][2].intersection(peak_range)
                
            #If there is an annotated site, add to new dictionary
            if len(five_prime) > 0:
                five_found_counter += 1
                five_prime = next(iter(five_prime))
                new_dict[gene][0].append(five_prime)
                new_dict[gene][2].append("5' splice site")
            elif len(three_prime) > 0:
                three_prime = next(iter(three_prime))
                new_dict[gene][0].append(three_prime)
                three_found_counter += 1
                new_dict[gene][2].append("3' splice site")
                
            #Otherwise, annotate as an unknown peak
            elif len(five_prime) == 0 and len(three_prime) == 0:
                new_counter += 1
                new_dict[gene][0].append(peak_tuples[n][0])
                new_dict[gene][2].append("Unknown")
                
        #Go through annotated splice sites and add those that were not found to the dictionary
        peaks = []
        for peak in new_peak_dict[gene][0]:
            peaks.append(peak[0])
            
        for five_site in new_peak_dict[gene][1]:
            site_range = range(five_site-3, five_site+4)
            match = set(site_range).intersection(peaks)
            if len(match) == 0:
                new_dict[gene][0].append(five_site)
                new_dict[gene][1].append(0)
                new_dict[gene][2].append("5' splice site")
                counter += 1
            
        for three_site in new_peak_dict[gene][2]:
            site_range = range(three_site-3, three_site+4)
            match = set(site_range).intersection(peaks)
            if len(match) == 0:
                new_dict[gene][0].append(three_site)
                new_dict[gene][1].append(0)
                new_dict[gene][2].append("3' splice site")
                counter += 1

        #Count all the splice sites
        for five_prime in new_peak_dict[gene][1]:
            five_counter += 1           
        for three_prime in new_peak_dict[gene][2]:
            three_counter += 1
            
    #Write out peaks with classifications to a tab delimited file
    transcript_counter = 0
    with open("{}_peak_picking.txt".format(prefix), "w") as fout:
        fout.write("Transcript\t Location\t Peak height\t P-value\t Classification\n")
        for gene, values in new_dict.iteritems():
            transcript_counter += 1
            n = 0
            while n < len(values[0]):
                line_list = [gene, values[0][n], values[1][n], values[2][n], "\n"]
                line_list = map(str, line_list)
                line = "\t".join(line_list)
                fout.write(line)
                n+=1
    
    #Print out the final counts for each peak type and return the new dictionary
    fivep_percent = float(five_found_counter)/five_counter*100
    threep_percent = float(three_found_counter)/three_counter*100
    
    print "\n"
    print prefix
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
    
def compare_peaks(dict_list, gff3_file, transcript_dict, all_transcripts=False, organism=None):
    transcript_set = set()
    splice_site_dict, intron_flag = list_splice_sites(gff3_file, organism=organism)
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
        if organism == 'pombe':
            fivep_total += len(splice_site_dict[transcript+".1"][0])
            threep_total += len(splice_site_dict[transcript+".1"][1])
            chrom = transcript_dict[transcript+".1"][3]
        else:
            fivep_total += len(splice_site_dict[transcript+"T0"][0])
            threep_total += len(splice_site_dict[transcript+"T0"][1])
            chrom = transcript_dict[transcript+"T0"][3]
        new_dict[transcript] = [[],[],[]]


        peak_list1 = dict_list[0][transcript][0]
        peak_list2 = dict_list[1][transcript][0]
        new_set1 = set((chrom, x) for x in peak_list1)
        new_set2 = set((chrom, x) for x in peak_list2)
        peak_set = new_set1.intersection(new_set2)

        for peak_dict in dict_list:
            i = 0
            for i in range(len(peak_dict[transcript][0])):
                peak = peak_dict[transcript][0][i]
                height = peak_dict[transcript][1][i]
                site_class = peak_dict[transcript][2][i]
                if peak in new_dict[transcript][0] or height == 0:
                    continue
                elif (chrom, peak) in peak_set:
                    new_dict[transcript][0].append(peak)
                    new_dict[transcript][1].append(height)
                    new_dict[transcript][2].append(site_class)
                    if site_class == "5' splice site":
                        fivep_found += 1
                    elif site_class == "3' splice site":
                        threep_found += 1
                    elif site_class == "Unknown":
                        new_found += 1
    
    fivep_percent = float(fivep_found)/fivep_total*100
    threep_percent = float(threep_found)/threep_total*100
    
    print "\n"
    #print dict_list
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
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N', 'Y':'R', 'R':'Y'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
        return complement(s[::-1])
    
def find_new_peak_sequence(fasta_file, gff3_file, peak_dict1, peak_dict2=None, prefix="Peak", organism=None):
    with open("{0}_sequences.txt".format(prefix), "w") as fout:
        line = "transcript\t chromosome\t peak type\t coordinate\t peak height\t strand\t sequence\t extended sequence\t looks like\n"
        fout.write(line)

    #Read fasta file for chromosome into list
    fasta_dict = {}
    letters = ['A','C','G','T','N']
    n = 0
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">") and line[1] not in ['M','A']:
                if organism == 'pombe':
                    chr_rom = line.split(' ')[0].strip('>')
                    rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
                    if chr_rom not in rom_lat: continue
                    chr_num = rom_lat[chr_rom]
                else:
                    chr_num = line[1:-1].strip()
                n = chr_num[3:]
                print n
                fasta_dict[chr_num] = str()
            elif line[0] in letters and "chr"+str(n) in fasta_dict:
                fasta_dict["chr"+str(n)] = fasta_dict["chr"+str(n)]+line.strip()

    fasta_dict = collections.OrderedDict(sorted(fasta_dict.items()))
    
    
    #Build transcript dictionary: keys are CNAG, values are [start, end, strand, chromosome]
    transcript_dict = build_transcript_dict(gff3_file, organism=organism)
    
    #Find positions of new peaks - first file
    index_dict1 = {}
    unknown_counter = 0
    counter5 = 0
    counter3 = 0
    unk5 = 0
    unk3 = 0
    unk = 0
    
    vol5 = 0
    vol3 = 0
    vol5_unk = 0
    vol3_unk = 0
    vol_unk = 0
    
    for CNAG, peaks in peak_dict1.iteritems():
        if organism == 'pombe':
            chromosome = transcript_dict[CNAG+".1"][3]
            strand = transcript_dict[CNAG+".1"][2]
        else:
            chromosome = transcript_dict[CNAG+"T0"][3]
            strand = transcript_dict[CNAG+"T0"][2]
        
        i = 0
        for i in range(len(peaks[0])):
            #Find sequence surrounding peak
            position = peaks[0][i]
            classification = peaks[2][i]
            if strand == "+":
                sequence = fasta_dict[chromosome][(position-2):(position+6)]
                ext_seq = fasta_dict[chromosome][(position-14):(position+18)]
                if organism == 'pombe': 
                    sequence = fasta_dict[chromosome][(position-3):(position+5)]
                    ext_seq = fasta_dict[chromosome][(position-15):(position+17)]

            elif strand == "-":
                sequence = fasta_dict[chromosome][(position-5):(position+3)]
                ext_seq = fasta_dict[chromosome][(position-17):(position+15)]
                if organism == 'pombe':
                    sequence = fasta_dict[chromosome][(position-6):(position+2)]
                    ext_seq = fasta_dict[chromosome][(position-18):(position+14)]

                sequence = reverse_complement(sequence)
                ext_seq = reverse_complement(ext_seq)
                
            #Classify sequences as splice sites
            if sequence[3:5] == "GT" or sequence[3:5] == "GC":
                site_class = "5'"
                GUpos = 3
                if classification == "5' splice site":
                    counter5 += 1
                    vol5 += peaks[1][i]
            elif sequence[1:3] == "AG":
                site_class = "3'"
                AGpos = 1
                if classification == "3' splice site":
                    counter3 += 1
                    vol3 += peaks[1][i]
                    
            #elif sequence[3:5] != "GT" and sequence[3:5] != "GC" and "GT" in sequence:
            #    site_class = "5'"
            #    GUpos = sequence.index('GT')
            #    if classification == "5' splice site":
            #        counter5 += 1
            #        vol5 += peaks[1][i]
            #elif sequence[3:5] != "GT" and sequence[3:5] != "GC" and "AG" in sequence:
            #    site_class = "3'"
            #    AGpos = sequence.index('AG')
            #    if classification == "3' splice site":
            #        counter3 += 1
            #        vol3 += peaks[1][i]
            else:
                site_class = "Unknown"
                unknown_counter += 1
            
            #if site_class == "5' offset":
            #    if transcript_dict[CNAG+"T0"][2] == "+":
            #        sequence = fasta_dict[chromosome][(position+GUpos-5):(position+GUpos+3)]
            #        ext_seq = fasta_dict[chromosome][(position+GUpos-17):(position+GUpos+15)]
            #    elif transcript_dict[CNAG+"T0"][2] == "-":
            #        sequence = fasta_dict[chromosome][(position-GUpos-2):(position-GUpos+6)]
            #        ext_seq = fasta_dict[chromosome][(position+GUpos-14):(position+GUpos+18)]
            #        sequence = reverse_complement(sequence)
            #        ext_seq = reverse_complement(ext_seq)
            #elif site_class == "3' offset":
            #    if transcript_dict[CNAG+"T0"][2] == "+":
            #        sequence = fasta_dict[chromosome][(position+AGpos-3):(position+AGpos+5)]
            #        ext_seq = fasta_dict[chromosome][(position+GUpos-15):(position+GUpos+17)]
            #    elif transcript_dict[CNAG+"T0"][2] == "-":
            #        sequence = fasta_dict[chromosome][(position-AGpos-4):(position-AGpos+4)]
            #        ext_seq = fasta_dict[chromosome][(position+GUpos-16):(position+GUpos+16)]
            #        sequence = reverse_complement(sequence)
            #        ext_seq = reverse_complement(ext_seq)
                    
            if classification == "Unknown" and "5'" in site_class: 
                unk5 += 1
                vol5_unk += peaks[1][i]
            elif classification == "Unknown" and "3'" in site_class: 
                unk3 += 1
                vol3_unk += peaks[1][i]
            elif classification == "Unknown" and site_class == "Unknown": 
                unk += 1
                vol_unk += peaks[1][i]
                
            #Build index dictionary: [transcript, chromosome, known site?, peak position, peak height, sequence (-4 to +4), classification
            if CNAG in index_dict1:      
                index_dict1[CNAG].append([CNAG, chromosome, peaks[2][i], position, peaks[1][i], strand, sequence, ext_seq, site_class])
            else:
                index_dict1[CNAG] = []
                index_dict1[CNAG].append([CNAG, chromosome, peaks[2][i], position, peaks[1][i], strand, sequence, ext_seq, site_class])
    
    
    #Make output file
    print "Output location: {0}_sequences.txt".format(prefix)
    
    print "5' splice sites: "+str(counter5)
    print "3' splice sites: "+str(counter3)
    print "Sequences with no splice sites: "+str(unknown_counter)
    
    print "Of unannotated splice sites:"
    print str(unk5)+" 5' splice sites, "+str(unk3)+" 3' splice sites, "+str(unk)+" other sites"
    print "Peak volumes:"
    print "5' annotated = "+str(vol5)
    print "3' annotated = "+str(vol3)
    print "New 5' = "+str(vol5_unk)
    print "New 3' = "+str(vol3_unk)
    print "Unknown sites = "+str(vol_unk)
    

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

def peak_density(peak_dict, gff3_file, by_intron=True, by_exon=False, transcript_list=None):
    transcript_dict = build_transcript_dict(gff3_file)
    if transcript_list is not None:
        transcript_dict = {key: transcript_dict[key] for key in transcript_list}
        peak_dict = {key: peak_dict[key] for key in transcript_list if key in peak_dict}
    
    splice_site_dict, intron_flag = list_splice_sites(gff3_file)
    
    peaks_by_intron = {}
    peaks_by_exon = {}
    for transcript, peaks in peak_dict.iteritems():
        chrom = transcript_dict[transcript][3]
        strand = transcript_dict[transcript][2]
        fivep_sites = splice_site_dict[transcript][0]
        threep_sites = splice_site_dict[transcript][1]
        i=0
        for i in range(len(peaks[0])):
            n=0
            for n in range(len(fivep_sites)):
                if strand == "+":
                    intron = transcript+"-"+str(n+1)
                    length = threep_sites[n]-fivep_sites[n]
                    if intron not in peaks_by_intron:
                        peaks_by_intron[intron]=[length,0]
                    if peaks[0][i] > fivep_sites[n] and peaks[0][i] < threep_sites[n] and peaks[3][i] == "Unknown":
                        peaks_by_intron[intron][1] += 1
                elif strand == "-":
                    intron = transcript+"-"+str(len(fivep_sites)-n)
                    length = fivep_sites[n]-threep_sites[n]
                    if intron not in peaks_by_intron:
                        peaks_by_intron[intron]=[length,0]
                    if peaks[0][i] > threep_sites[n] and peaks[0][i] < fivep_sites[n] and peaks[3][i] == "Unknown":
                        peaks_by_intron[intron][1] += 1
    
            n=0
            for n in range(len(fivep_sites)+1):    
                if strand == "+":
                    exon = transcript+"-"+str(n)
                    if n == 0:
                        exon_length = fivep_sites[n]-transcript_dict[transcript][0]
                    elif n > 0 and n < len(fivep_sites):
                        exon_length = fivep_sites[n]-threep_sites[n-1]
                    elif n == len(fivep_sites):
                        exon_length = transcript_dict[transcript][1]-threep_sites[n-1]
                    if exon not in peaks_by_exon:
                        peaks_by_exon[exon]=[exon_length,0]
                    if peaks[3][i] == "Unknown":
                        if n == 0 and peaks[0][i] > transcript_dict[transcript][0] and peaks[0][i] < fivep_sites[n]:
                            peaks_by_exon[exon][1] += 1
                        elif n>0 and n<len(fivep_sites) and peaks[0][i]>threep_sites[n-1] and peaks[0][i]<fivep_sites[n]:
                            peaks_by_exon[exon][1] += 1
                        elif n == len(fivep_sites) and peaks[0][i]>threep_sites[n-1] and peaks[0][i]<transcript_dict[transcript][1]:
                            peaks_by_exon[exon][1] += 1
                if strand == "-":
                    exon = transcript+"-"+str(len(fivep_sites)-n)
                    if n == 0:
                        exon_length = threep_sites[n]-transcript_dict[transcript][0]
                    elif n>0 and n<len(fivep_sites):
                        exon_length = threep_sites[n]-fivep_sites[n-1]
                    elif n == len(fivep_sites):
                        exon_length = transcript_dict[transcript][1]-fivep_sites[n-1]
                    if exon not in peaks_by_exon:
                        peaks_by_exon[exon]=[exon_length,0]
                    if peaks[3][i] == "Unknown":
                        if n == 0 and peaks[0][i] > transcript_dict[transcript][0] and peaks[0][i] < threep_sites[n]:
                            peaks_by_exon[exon][1] += 1
                        elif n>0 and n<len(fivep_sites) and peaks[0][i]>fivep_sites[n-1] and peaks[0][i]<threep_sites[n]:
                            peaks_by_exon[exon][1] += 1
                        elif n == len(fivep_sites) and peaks[0][i]>fivep_sites[n-1] and peaks[0][i] <transcript_dict[transcript][1]:
                            peaks_by_exon[exon][1] += 1

        
    
    print len(peaks_by_intron)
    print len(peaks_by_exon)
    lengths = []
    num_peaks = []
    if by_intron==True:
        for intron, values in peaks_by_intron.iteritems():
            lengths.append(values[0])
            num_peaks.append(values[1])
    elif by_exon==True:
        for exon, values in peaks_by_exon.iteritems():
            lengths.append(values[0])
            num_peaks.append(values[1])
        
    fig1 = pyplot.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(lengths, num_peaks, c='royalblue', alpha = 0.5)
    pyplot.show()
    return fig1

###################################################################################################
## Output a dictionary of known splice sites                 ##
## Takes peak dictionary from compare_peaks and transcript dictionary from build_transcript_dict ##
###################################################################################################

def make_gff3_dict(gff3_file):
    gff3_dict = {}
    with open(gff3_file, "r") as fin:
        for line in fin:
            if line.startswith('chr'):
                line_list = line.split("\t")
                if line_list[2] == "mRNA":
                    gene = line_list[8].split("=")[1].split("T")[0]
                    chrom = line_list[0]
                    start = int(line_list[3])
                    stop = int(line_list[4])
                    strand = line_list[6]
                    gff3_dict[gene] = [chrom, start, stop, strand, []]
                elif line_list[2] == "exon":
                    gene = line_list[8].split("=")[2].split("T")[0]
                    start = int(line_list[3])
                    stop = int(line_list[4])
                    if start != gff3_dict[gene][1] and start != gff3_dict[gene][2]:
                        gff3_dict[gene][4].append(start)
                    if stop != gff3_dict[gene][1] and stop != gff3_dict[gene][2]:
                        gff3_dict[gene][4].append(stop)
    for gene, values in gff3_dict.iteritems():
        if values[1] in values[4]: 
            print "Start"
            print gene
            print values[1]
            print strand
        elif values[2] in values[4]: 
            print "Stop"
            print gene
            print values[2]
            print strand
    return gff3_dict

###################################################################################################
## Output putative spliceosomal cleavages to new gff3 file for use with crunchBAM                ##
## Takes peak dictionary from compare_peaks and transcript dictionary from build_transcript_dict ##
###################################################################################################

def new_peaks_gff3(peak_dict, transcript_dict, prefix='NEW_PEAKS', organism=None):
    gff3_dict = {}
    with open("{0}.gff3".format(prefix), "w") as fout:
        fout.write("")
    
    alias_base = 7e16 + random.randint(0,1e11)
    alias_base = int(alias_base)
    n = 0
    for gene, peaks in peak_dict.iteritems():
        alias_base += 1
        peak_list = []
        if organism == 'pombe':
            chrom = transcript_dict[gene+'.1'][3]
            start = transcript_dict[gene+'.1'][0]
            stop = transcript_dict[gene+'.1'][1]
            strand = transcript_dict[gene+'.1'][2]
            CDS_starts = transcript_dict[gene+'.1'][4]
            CDS_stops = transcript_dict[gene+'.1'][5]
        else:
            chrom = transcript_dict[gene+'T0'][3]
            start = transcript_dict[gene+'T0'][0]
            stop = transcript_dict[gene+'T0'][1]
            strand = transcript_dict[gene+'T0'][2]
            CDS_starts = transcript_dict[gene+'T0'][4]
            CDS_stops = transcript_dict[gene+'T0'][5]
        gff3_dict[gene] = [chrom, start, stop, strand, []]

        n = 0
        for n in range(len(peaks[0])):
            if peaks[2][n] == "Unknown":
                gff3_dict[gene][4].append(peaks[0][n])
                peak_list.append(peaks[0][n])
                
        if len(peak_list) > 0:
            start = str(start)
            stop = str(stop)
            with open("{0}.gff3".format(prefix), "a") as fout:
                gene_line_list = ["\n"+chrom, prefix, "gene", start, stop, ".", strand.strip(), ".", "ID="+gene+";Alias="+('%d' % alias_base)+";Name=Unknown\n"]
                gene_line = "\t".join(gene_line_list)
                fout.write(gene_line)
            
                n = 0
                for n in range(len(peak_list)):
                    alias_base += 1
                    mRNA_line_list = [chrom, prefix, "mRNA", start, stop, ".", strand.strip(), ".", "ID="+gene+"T"+str(n)+";Alias="+('%d' % alias_base)+";Parent="+gene+";Name=Unknown\n"]
                    mRNA_line = "\t".join(mRNA_line_list)
                    fout.write(mRNA_line)
            
                    alias_base += 1
                    if len(CDS_starts) > 0:
                        if strand == '+':
                            CDS_line_list = [chrom, prefix, "CDS", CDS_starts[0], CDS_stops[-1], ".", strand.strip(), ".", "ID="+gene+"T0.cds;Parent="+gene+"T"+str(n)+"\n"]
                        elif strand == '-':
                            CDS_line_list = [chrom, prefix, "CDS", CDS_starts[-1], CDS_stops[0],".", strand.strip(), ".", "ID="+gene+"T"+str(n)+".cds;Parent="+gene+"T"+str(n)+"\n"]
                        CDS_line = "\t".join(CDS_line_list)
                        fout.write(CDS_line)
            
                    alias_base += (1+n)
                    if strand == "+":
                        peak_line_list = [chrom, prefix, "exon", start, str(peak_list[n]), ".", strand.strip(), ".", "ID="+('%d' % alias_base)+";Parent="+gene+"T"+str(n)+"\n"]
                    elif strand == "-":
                        peak_line_list = [chrom, prefix, "exon", str(peak_list[n]), stop, ".", strand.strip(), ".", "ID="+('%d' % alias_base)+";Parent="+gene+"T"+str(n)+"\n"]
                    if int(peak_line_list[4])-int(peak_line_list[3]) <= 0: print gene+" "+strand
                
                    peak_line = "\t".join(peak_line_list)
                    fout.write(peak_line)
    return gff3_dict

def count_reads_at_splice_sites(gff3_dict, transcript_dict, bedgraph_list, sort_bedgraph=True):
    bedgraph_dict_list = []
    if sort_bedgraph == True:
        for bedgraph in bedgraph_list:
            bedgraph_dict_list.append(build_bedgraph_dict(transcript_dict, bedgraph))
    elif sort_bedgraph == False:
        for bedgraph in bedgraph_list:
            bedgraph_dict_list.append(read_CNAGsort_bedgraph(bedgraph))
            
    index_tuples = []
    for gene, coords in gff3_dict.iteritems():
        n = 0
        for n in range(len(coords[4])):
            index_tuples.append((gene+'T0', coords[4][n]))
    
    index = pd.MultiIndex.from_tuples(index_tuples)
    columns = [('Transcript0','Transcript0'),('Exon0','Exon0')]
    for bg in bedgraph_list:
        columns.append((bg.split("/")[-1].split("_")[0],'5prime'))
    for bg in bedgraph_list:
        columns.append((bg.split("/")[-1].split("_")[0],'3prime'))
    columns = pd.MultiIndex.from_tuples(columns)
    count_df = pd.DataFrame(index=index, columns=columns)
    
    n=0
    for n in range(len(bedgraph_dict_list)):
        bedgraph_dict = bedgraph_dict_list[n]
        for gene, values in bedgraph_dict.iteritems():
            if gene.split('T')[0] in gff3_dict:
                peak_set = set(values[0]).intersection(gff3_dict[gene.split('T')[0]][4])
                no_bg_set = set(gff3_dict[gene.split('T')[0]][4]).difference(values[0])
                peak_set.update(no_bg_set)
                for peak in peak_set:
                    peak_size = 0
                    peak_range = range(peak-2, peak+2)
                    for coord in peak_range:
                        if coord in values[0]:
                            a = values[0].index(coord)
                            peak_size += values[1][a]
                    count_df.loc[(gene.split('T')[0]+('T0'), peak),columns[0]] = gene.split('T')[0]+('T0')
                    count_df.loc[(gene.split('T')[0]+('T0'), peak),columns[1]] = peak
                    count_df.loc[(gene.split('T')[0]+('T0'), peak),columns[n+2]] = peak_size
                    
    return count_df

def peak_pipeline(gff3_file, fasta_file, piranha_bedgraph_tuples, conditions=1, prefix='', organism=None):
    #Build dictionary of peaks picked by piranha. Requires 5' limited bedgraph file which will be sorted by gene and saved
    transcript_dict = build_transcript_dict(gff3_file, organism=organism)
    peak_dicts = []
    bedgraph_list = []
    for piranha_output, bedgraph in piranha_bedgraph_tuples:
        bedgraph_list.append(bedgraph)
        peak_dicts.append(build_piranha_dict(transcript_dict, piranha_output, sort_bedgraph=False, sorted_bedgraph=bedgraph, max_p_value=1))
        
    #Classify peaks based on whether they occur at annotated splice sites
    new_peaks = []
    for peak_dict in peak_dicts:
        new_peaks.append(check_splice_sites(peak_dict, gff3_file, prefix=prefix, organism=organism))
        
    #Compare peaks between samples. Report only those occuring in both samples
    both_dict_list = []
    if conditions == 1:
        all_dict = compare_peaks(new_peaks, gff3_file, transcript_dict, all_transcripts=False, organism=organism)
    else:
        print str(conditions)+" conditions"
        n = 0
        while n < conditions:
            both_dict_list.append(compare_peaks([new_peaks[n], new_peaks[n+1]], gff3_file, transcript_dict, all_transcripts=False, organism=organism))
            n += 2
        all_dict = both_dict_list[0].copy()
        a = 1
        for a in range(len(both_dict_list)):
            all_dict.update(both_dict_list[a])
    
    #Get sequences surrounding all peaks and classify based on whether they contain a splice site sequence
    find_new_peak_sequence(fasta_file, gff3_file, all_dict, prefix=prefix, organism=organism)
    
    #Create a dictionary from gff3 file based on discovered peaks
    gff3_dict = new_peaks_gff3(all_dict, transcript_dict, prefix=prefix, organism=organism)
    
    #Create the final dataframe and save as a tsv. This can now be used with the SPTables module
    new_peaks_df = count_reads_at_splice_sites(gff3_dict, transcript_dict, bedgraph_list, sort_bedgraph=False)
    new_peaks_df.to_csv('{0}_called_peaks.tsv'.format(prefix), sep='\t')
    
    return new_peaks_df

def annotated_pipeline(gff3_file, sorted_bedgraph_list, prefix='', count_reads=True, count_file=None):
    transcript_dict = build_transcript_dict(gff3_file)
    gff3_dict = make_gff3_dict(gff3_file)   
    if count_reads is True:
        df = count_reads_at_splice_sites(gff3_dict, transcript_dict, sorted_bedgraph_list, sort_bedgraph=False)
        df.to_csv('{0}_reads_at_ss.tsv'.format(prefix), sep='\t')
    elif count_reads is False:
        df = SPTables.build_tables(count_file, header=[0,1], skiprows=[2],  multiIndex=True)
        df.drop(('Transcript1','Transcript1'), axis=1, inplace=True)
        df.drop(('Exon2','Exon2'), axis=1, inplace=True)
    return df
    

def normalize_and_plot_peaks(peak_df, gobs_count_reads_in_transcripts, transcript_len_file, control_file, rep_tuples, mutant='condition 2', cutoff = None):
    #transcript_len_file for crypto: "/home/jordan/CodeBase/RNA-is-awesome/GENOMES/H99_transcript_lengths.txt"
    if type(gobs_count_reads_in_transcripts) == str:
        totals_df = SPTables.build_tables(gobs_count_reads_in_transcripts)
    elif type(gobs_count_reads_in_transcripts) == list:
        totals_df = SPTables.build_tables(gobs_count_reads_in_transcripts[0], file2=gobs_count_reads_in_transcripts[1])
    transcript_lengths = SPTables.build_tables(transcript_len_file)
    controls = SPTables.build_tables(control_file, int_index=True)
    df = SPTables.normalize_AtoB(peak_df, totals_df, transcript_lengths, controls)
    index_tuples = zip(df[('Transcript0','Transcript0')].tolist(), df[('Exon0','Exon0')].tolist())
    df.index = pd.MultiIndex.from_tuples(index_tuples)
    
    #Average replicates if indicated
    new_df = pd.DataFrame(index=df.index)
    for n in range(len(rep_tuples)):
        new_df[(rep_tuples[n][0],'5prime Normalized')] = df[(rep_tuples[n][0],'5prime Normalized')]
        new_df[(rep_tuples[n][1],'5prime Normalized')] = df[(rep_tuples[n][1],'5prime Normalized')]
        new_df[(rep_tuples[n][0]+'-A','5prime')] = df[(rep_tuples[n][0]+'-A','5prime')]
        new_df[(rep_tuples[n][0]+'-B','Total')] = df[(rep_tuples[n][0]+'-B','Total')]
        new_df[(rep_tuples[n][1]+'-A','5prime')] = df[(rep_tuples[n][1]+'-A','5prime')]
        new_df[(rep_tuples[n][1]+'-B','Total')] = df[(rep_tuples[n][1]+'-B','Total')]
        new_df[('avg'+str(n+1),'5prime Normalized')] = df[[(rep_tuples[n][0],'5prime Normalized'), (rep_tuples[n][1],'5prime Normalized')]].mean(axis=1)
        new_df[('A-avg'+str(n+1),'5prime')] = df[[(rep_tuples[n][0]+'-A','5prime'),(rep_tuples[n][1]+'-A','5prime')]].mean(axis=1)
        new_df[('B-avg'+str(n+1),'Total')] = df[[(rep_tuples[n][0]+'-B','Total'),(rep_tuples[n][1]+'-B','Total')]].mean(axis=1)
         
    if cutoff is not None and len(rep_tuples) > 1:
        new_df = new_df[new_df[('B-avg2','Total')] > cutoff]
    #elif cutoff is not None and len(rep_tuples) == 1:

        
    #Sort peaks by whether they increase or decrease between conditions
    columns = pd.MultiIndex.from_tuples(new_df.columns)
    increase_df = pd.DataFrame(index=new_df.index, columns=columns)
    decrease_df = pd.DataFrame(index=new_df.index, columns=columns)
    other_df = pd.DataFrame(index=new_df.index, columns=columns)
    
    if len(rep_tuples) > 1:            
        increase_df = new_df[new_df[(rep_tuples[1][0],'5prime Normalized')] > new_df[(rep_tuples[0][0],'5prime Normalized')]*3]
        increase_df = increase_df[increase_df[(rep_tuples[1][1],'5prime Normalized')] > increase_df[(rep_tuples[0][1],'5prime Normalized')]*3]
        decrease_df = new_df[new_df[(rep_tuples[1][0],'5prime Normalized')] < new_df[(rep_tuples[0][0],'5prime Normalized')]*0.33]
        decrease_df = decrease_df[decrease_df[(rep_tuples[1][1],'5prime Normalized')] < decrease_df[(rep_tuples[0][1],'5prime Normalized')]*0.33]
        other_df = new_df[new_df[(rep_tuples[1][0],'5prime Normalized')] <= new_df[(rep_tuples[0][0],'5prime Normalized')]*3]
        other_df = other_df[other_df[(rep_tuples[1][1],'5prime Normalized')] <= other_df[(rep_tuples[0][1],'5prime Normalized')]*3]
        other_df = other_df[other_df[(rep_tuples[1][0],'5prime Normalized')] >= other_df[(rep_tuples[0][0],'5prime Normalized')]*0.33]
        other_df = other_df[other_df[(rep_tuples[1][1],'5prime Normalized')] >= other_df[(rep_tuples[0][1],'5prime Normalized')]*0.33]

    #increase_df = increase_df[increase_df[('avg2','5prime Normalized')].map(str) != 'nan']
    #decrease_df = decrease_df[decrease_df[('avg2','5prime Normalized')].map(str) != 'nan']
    #other_df = other_df[other_df[('avg2','5prime Normalized')].map(str) != 'nan']
        
    #Make lists for plotting

        other1 = other_df[('avg1','5prime Normalized')].tolist()
        other2 = other_df[('avg2','5prime Normalized')].tolist()
        inc1 = increase_df[('avg1','5prime Normalized')].tolist()
        inc2 = increase_df[('avg2','5prime Normalized')].tolist()
        dec1 = decrease_df[('avg1','5prime Normalized')].tolist()
        dec2 = decrease_df[('avg2','5prime Normalized')].tolist()
        
    elif len(rep_tuples) == 1:
        increase_df = new_df[new_df[(rep_tuples[0][1],'5prime Normalized')] > new_df[(rep_tuples[0][0],'5prime Normalized')]*3]
        decrease_df = new_df[new_df[(rep_tuples[0][1],'5prime Normalized')] < new_df[(rep_tuples[0][0],'5prime Normalized')]*0.33]
        other_df = new_df[new_df[(rep_tuples[0][1],'5prime Normalized')] <= new_df[(rep_tuples[0][0],'5prime Normalized')]*3]
        other_df = other_df[other_df[(rep_tuples[0][1],'5prime Normalized')] >= other_df[(rep_tuples[0][0],'5prime Normalized')]*0.33]
        
        other1 = other_df[(rep_tuples[0][0],'5prime Normalized')].tolist()
        other2 = other_df[(rep_tuples[0][1],'5prime Normalized')].tolist()
        inc1 = increase_df[(rep_tuples[0][0],'5prime Normalized')].tolist()
        inc2 = increase_df[(rep_tuples[0][1],'5prime Normalized')].tolist()
        dec1 = decrease_df[(rep_tuples[0][0],'5prime Normalized')].tolist()
        dec2 = decrease_df[(rep_tuples[0][1],'5prime Normalized')].tolist()

    else: print "Please provide the names of replicates in the rep_tuples argument"
    
    #Convert lists to log values for plotting
    other1 = convert_to_log(other1)
    other2 = convert_to_log(other2)
    inc1 = convert_to_log(inc1)
    inc2 = convert_to_log(inc2)
    dec1 = convert_to_log(dec1)
    dec2 = convert_to_log(dec2)

    
    #Add up some lists for statistical tests
    all_log = other1+inc1+dec1+other2+inc2+dec2
    all1 = other1+inc1+dec1
    all2 = other2+inc2+dec2
    print "\n2 sample T-test"
    t, p = stats.ttest_ind(numpy.array(all2), numpy.array(all1))
    print "T-value "+str(t)
    print "P-value: "+str(p)+"\n"
    
    #Scatter plot
    scatter = pyplot.figure()
    ax1 = scatter.add_subplot(111)
    ax1.scatter(other1, other2, c='0.5', label='No change', alpha = 0.5, edgecolor='0.3')
    ax1.scatter(inc1, inc2, c='coral', label='Increased in '+mutant, alpha = 0.5, edgecolor='coral')
    ax1.scatter(dec1, dec2, c='royalblue', label='Decreased in '+mutant, alpha = 0.5, edgecolor='darkslateblue')
    xmax = numpy.nanmax(all_log)+0.5
    xmin = numpy.nanmin(all_log)-0.5
    ax1.set_ylim([xmin,xmax])
    ax1.set_xlim([xmin,xmax])
    ymax = ax1.get_ylim()
    ax1.legend(loc=4)
    ax1.plot([xmin, xmax], [xmin, xmax], ls="--", c=".3", lw=1.5)
    
    print "Linear regression:"
    slope, intercept, r_value, p_value, std_err = stats.linregress(all1,all2)
    print 'Slope: '+str(slope)
    print 'Intercept: '+str(intercept)
    print 'R^2 value: '+str(r_value**2)
    print 'P value: '+str(p_value)
    print 'Standard error '+str(std_err)
    
    return (increase_df, scatter)
        
def convert_to_log(num_list):
    log_list = [numpy.NaN if x == 0 else log(x, 10) for x in num_list]
    return log_list

def analyze_peaks(peak_tsv_file, gff3_file):
    df = pd.read_csv(peak_tsv_file, sep='\t')
    print df.columns
    transcripts = list(set(df['transcript'].tolist()))
    print len(transcripts)
    transcripts_T0 = [x+'T0' for x in transcripts]
    transcript_dict = build_transcript_dict(gff3_file)
    splice_site_dict, intron_flag = list_splice_sites(gff3_file, gene_list = transcripts_T0)
    print len(splice_site_dict)
    
    index = pd.MultiIndex(levels=[[],[]], labels=[[],[]], names=[u'transcript', u'interval'])
    new_df = pd.DataFrame(columns=['chromosome','# novel peaks', 'peak coords', 'peak heights', 'annotated sites', 'interval length'], index=index)
    
    for transcript in transcripts:
        #print transcript
        tx_df = df[df['transcript'] == transcript]
        tx_df = tx_df.reset_index()
        tx_df['intron/exon'] = 'Unknown'
        intervals = []
        intervals.append(transcript_dict[transcript+'T0'][0])
        intervals.append(transcript_dict[transcript+'T0'][1])
        intervals = intervals + splice_site_dict[transcript+'T0'][0]
        intervals = intervals + splice_site_dict[transcript+'T0'][1]
        intervals.sort()

        a = 0
        for a in range(len(intervals)-1):
            interval = a
            interval_length = intervals[a+1] - intervals[a]
            if transcript_dict[transcript+'T0'][2] == '-':
                interval = len(intervals)-interval
                
            peak_count = 0
            peak_coords = []
            peak_heights = []
            annotated_sites = []
            n=0
            for n in range(len(tx_df)): 
                if tx_df[' coordinate'][n] < intervals[a+1] and tx_df[' coordinate'][n] > intervals[a]:
                    peak_count += 1
                    peak_coords.append(tx_df[' coordinate'][n])
                    peak_heights.append(tx_df[' peak height'][n])
                elif tx_df[' coordinate'][n] == intervals[a+1] or tx_df[' coordinate'][n] == intervals[a]:
                    annotated_sites.append(tx_df[' coordinate'][n])

            if interval % 2 == 0:
                interval = 'exon '+str(interval/2)
            else:
                interval = 'intron '+str(interval/2+1)
            if peak_count != 0:
                new_df.loc[(transcript,interval),] = [transcript_dict[transcript+'T0'][3], peak_count, ','.join(str(x) for x in peak_coords), ','.join(str(x) for x in peak_heights), ','.join(str(x) for x in annotated_sites), interval_length]

            #new_df.loc[(transcript, interval), 'chromosome'] = transcript_dict[transcript+'T0'][3]
            #new_df.loc[(transcript, interval), '# novel peaks'] = peak_count
            #new_df.loc[(transcript, interval), 'peak coords'] = ','.join(str(x) for x in peak_coords)
            #new_df.loc[(transcript, interval), 'peak heights'] = ','.join(str(x) for x in peak_heights)
            #new_df.loc[(transcript, interval), 'annotated sites'] = ','.join(str(x) for x in annotated_sites)
            
    new_df.to_csv('{0}_by_interval.tsv'.format(peak_tsv_file.split('.')[-2]), sep='\t')
    return new_df
            
def check_annotation(peak_pipeline_output, gff3_file, organism=None):
    transcript_dict = build_transcript_dict(gff3_file, organism=organism)
    splice_site_dict, intron_flag = list_splice_sites(gff3_file, organism=organism)
    peak_df = pd.read_csv(peak_pipeline_output, sep='\t')
    
    misannotated = {}
    for transcript, sites in splice_site_dict.iteritems():
        misannotated[transcript] = []
        all_sites = sites[0]+sites[1]
        temp_df = peak_df[peak_df['transcript'] == transcript[:-2]]
        #print temp_df.columns
        peak_list = zip(temp_df[' coordinate'].tolist(), temp_df[' looks like'].tolist())
        distance_dict = {}
        if len(peak_list) > 0:
            for peak in peak_list:
                peak_coord = int(peak[0])
                peak_seq = peak[1]
                if peak_coord in all_sites: continue
                if peak_seq == "5'":
                    for site in sites[0]:
                        distance_dict[site] = abs(peak_coord-site)
                elif peak_seq == "3'":
                    for site in sites[1]:
                        distance_dict[site] = abs(peak_coord-site)
                sorted_dd = sorted(distance_dict.items(), key=operator.itemgetter(1))
                try:
                    if sorted_dd[0][0] in peak_list: continue
                    misannotated[transcript].append([peak, sorted_dd[0][0], sorted_dd[0][1]])
                except IndexError:
                    print transcript
                    print peak
    
    misannotated = dict((k, v) for k, v in misannotated.iteritems() if v)
    return misannotated
            
    
#####################################################################################################
## New simplified code for processing peaks from Jessica Li                                        ##
#####################################################################################################

#Function to read peak output file
def CP_peaks_by_gene(fin, transcript_dict, cutoff=5):
    genes_by_chr = {}
    for tx, info in transcript_dict.iteritems():
        if info[3] not in genes_by_chr:
            genes_by_chr[info[3]] = []
        genes_by_chr[info[3]].append(tx)
    
    rom_lat = {'I':'chr1','II':'chr2','III':'chr3'}
    
    peak_count = 0
    line_count = 0
    peaks_by_gene = {}
    strand_dict = {'0':'-','1':'+'}
    with open(fin,'r') as f:
        for line in f:
            data = line.split('\t')
            if len(data[0]) > 0 and len(data) > 3:
                chrom = data[0]
                peak = int(data[1])
                strand = strand_dict[data[2]]
                peak_height = float(data[3])
            elif len(data[0]) > 0 and len(data) == 3:
                chrom = data[0]
                if '_' in chrom:
                    chrom = 'chr'+chrom.split('_')[-1]
                peak = int(data[1])
                peak_height = float(data[2])
                strand = None
            else:
                chrom = data[1]
                peak = int(data[2])
                strand = strand_dict[data[3]]
                peak_height = float(data[4])

            if chrom in rom_lat: chrom = rom_lat[chrom]
            
            if peak_height >= cutoff:
                line_count += 1
                if chrom in genes_by_chr:
                    tx_list = genes_by_chr[chrom]
                    for tx in tx_list:
                        start = transcript_dict[tx][0]
                        end = transcript_dict[tx][1]
                        if peak > start and peak < end and (strand == transcript_dict[tx][2] or strand is None):
                            peak_count += 1
                            if tx not in peaks_by_gene:
                                peaks_by_gene[tx] = []
                            peaks_by_gene[tx].append([peak,peak_height,strand])
                            
                            ## Bit of code to check for a neighboring peak and remove it if there is a neighbor.
                            #else:
                            #    n=0
                            #    neighbor=False
                            #    for n in range(len(peaks_by_gene[tx])):
                            #        if abs(peak-peaks_by_gene[tx][n][0]) <= 2:
                            #            neighbor=True
                            #            if peak_height >= peaks_by_gene[tx][n][1]:
                            #                del peaks_by_gene[tx][n]
                            #                peaks_by_gene[tx].append([peak,peak_height,strand])
                            #            elif peak_height < peaks_by_gene[tx][n]:
                            #                pass
                            #    if neighbor == False:
                            #        peaks_by_gene[tx].append([peak,peak_height,strand])
    print peak_count
    return peaks_by_gene
 
#Function to compare untagged and 2 replicates and pick peaks that are in both but not in untagged
def CP_compare_reps(untagged, tagged1, tagged2):
    #print len(tagged1)
    #print len(tagged2)
    tx_list = list(set(tagged1.keys()).intersection(tagged2.keys()))
    #print len(tx_list)
    new_peak_dict = {}
    peak_count = 0
    for tx in tx_list:
        new_peak_dict[tx] = []
        for peak,peak_height,strand in tagged1[tx]:
            #if peak_height > 0.05*max(zip(*tagged1[tx])[1]):
            if peak in zip(*tagged2[tx])[0]:
                if tx not in untagged or peak not in zip(*untagged[tx])[0]:
                    new_peak_dict[tx].append([peak,peak_height,strand])
                    peak_count += 1
    print peak_count
    return new_peak_dict

#Function to check peaks against annotation
def CP_compare_to_annotation(peaks, ss_dict, transcript_dict):
    five_count = 0
    three_count = 0
    other_count = 0
    intronic_count = 0
    peak_count = 0
    for tx, peak_list in peaks.iteritems():
        peak_count += len(peak_list)   
    print peak_count
    compare_df = pd.DataFrame(index = range(peak_count+1), columns=['transcript','chromosome','strand','position','height','type'])
    n=0
    for tx, info in ss_dict.iteritems():
        chrom = transcript_dict[tx][3]
        if tx in peaks:
            if len(info[0]) > 0 and len(peaks[tx]) > 0:
                for peak, height, strand in peaks[tx]:
                    if strand is None:
                        strand = transcript_dict[tx][2]
                    peak_range = range(peak-3,peak+3)
                    try:
                        annotated = False
                        for pos in peak_range:
                        #if peak-1 in info[0] or peak in info[0] or peak+1 in info[0]:
                            if pos in info[0]:
                                five_count += 1
                                compare_df.ix[n] = [tx[:-2], chrom, strand, peak, height, "5prime"]
                                annotated = True
                                break
                            #print [tx[:-2], chrom, strand, peak, height, "5prime"]
                        #elif peak-1 in info[1] or peak in info[1] or peak+1 in info[1]:
                            elif pos in info[1]:
                                three_count += 1
                                compare_df.ix[n] = [tx[:-2], chrom, strand, peak, height, "3prime"]
                                annotated = True
                                break
                        if annotated is False:
                            other_count += 1
                            intron_flag = False
                            m=0
                            for m in range(len(info[0])):
                                if peak > info[0][m] and peak < info[1][m]:
                                    intronic_count += 1
                                    intron_flag = True
                                    break
                            if intron_flag is True:
                                compare_df.ix[n] = [tx[:-2], chrom, strand, peak, height, "intronic"]
                            else:
                                compare_df.ix[n] = [tx[:-2], chrom, strand, peak, height, "other"]
                    except IndexError:
                        print tx
                        print n
                        print peak_count
                        print len(info)
                    n+=1
    print "5prime annotated sites: "+str(five_count)
    print "3prime annotated sites: "+str(three_count)
    print "Unpredicted peaks: "+str(other_count)
    print "Unpredicted peaks in introns: "+str(intronic_count)
    compare_df.dropna(how='all',inplace=True)
    return compare_df


def collapse_unpredicted_peaks(df):
    tx_list = list(set(df['transcript'].tolist()))
    for tx in tx_list:
        tx_df = df[df['transcript'] == tx]
        tx_df = tx_df[tx_df['type'].isin(['other','intronic'])]
        index=tx_df.index
        n=0
        for n in range(len(index)):
            if tx_df['height'][index[n]] < 10:
                if index[n] in df.index:
                    df.drop(index[n], inplace=True)
            else:
                m=0
                for m in range(len(index)):
                    spacing = abs(tx_df['position'][index[n]]-tx_df['position'][index[m]])
                    if spacing > 0 and spacing <= 2:
                        if tx_df['height'][index[n]] > tx_df['height'][index[m]]:
                            if index[m] in df.index:
                                df.drop(index[m], inplace=True)
                        if tx_df['height'][index[n]] < tx_df['height'][index[m]]:
                            if index[n] in df.index:
                                df.drop(index[n], inplace=True)
                            break
                        else:
                            continue
    print "Number of unpredicted peaks after condensing:"
    print len(df[df['type'].isin(['other','intronic'])])
    print "Number of intronic peaks after condensing:"
    print len(df[df['type'] == 'intronic'])
    df.reset_index(inplace=True)
    return df

#Add sequences and check whether they're splice sites
def add_sequence_to_df(df, fa_dict):
    seq_df = df
    sequence = []
    looks_like = []
    for index, row in df.iterrows():
        if row['strand'] == '+':
            seq = fa_dict[row['chromosome']][row['position']-6:row['position']+6]
        elif row['strand'] == '-':
            seq = fa_dict[row['chromosome']][row['position']-7:row['position']+5]
            seq = reverse_complement(seq)
        sequence.append(seq)
        
        if row['type'] == '5prime':
            looks_like.append('5prime')
        elif row['type'] == '3prime':
            looks_like.append('3prime')
        else:
            if seq[6:8] == 'GT' or seq[6:8] == 'GC' or seq[6:8] == 'AT':
                looks_like.append(seq[6:8])
            elif seq[4:6] == 'AG' or seq[4:6] == 'AC':
                looks_like.append(seq[4:6])
            else:
                looks_like.append('')
    seq_df['sequence'] = sequence
    seq_df['looks like'] = looks_like
    return seq_df

def peak_to_seq_pipeline(untagged_peak_file, tagged1_peak_file, tagged2_peak_file, gff3, fasta, organism=None, cutoff=5):
    test_tx = 'SPCC16A11.10c'
    
    transcript_dict = build_transcript_dict(gff3, organism=organism)
    print "Finding peaks in transcripts..."
    print untagged_peak_file
    untagged = CP_peaks_by_gene(untagged_peak_file, transcript_dict, cutoff=cutoff)
    if test_tx in untagged: print untagged[test_tx]
    print tagged1_peak_file
    tagged1 = CP_peaks_by_gene(tagged1_peak_file, transcript_dict, cutoff=cutoff)
    if test_tx in tagged1: print untagged[test_tx]
    print tagged2_peak_file
    tagged2 = CP_peaks_by_gene(tagged2_peak_file, transcript_dict, cutoff=cutoff)
    if test_tx in tagged1: print untagged[test_tx]
    
    print "Comparing peaks between replicates..."
    peaks = CP_compare_reps(untagged, tagged1, tagged2)
    
    print "Checking peaks against annotation..."
    ss_dict, flag = list_splice_sites(gff3, organism=organism)
    peak_df = CP_compare_to_annotation(peaks, ss_dict, transcript_dict)
    peak_df = collapse_unpredicted_peaks(peak_df)
    
    if type(fasta) == str:
        fasta = SPScores.make_fasta_dict(fasta)
    print "Adding sequences..."
    peak_seq_df = add_sequence_to_df(peak_df, fasta)
    print "Completed"
    return peak_seq_df

def count_peak_types(df):
    other = df[df['type'].isin(['other','intronic'])]
    print len(other)
    other = other[other['looks like'] == '']
    print 'GT:'
    print len(df[df['looks like'] == 'GT'])
    print len(other[other['sequence'].str[5:7].str.contains('GT')])
    other = other[~other['sequence'].str[5:7].str.contains('GT')]
    print 'GC:'
    print len(df[df['looks like'] == 'GC'])
    print len(other[other['sequence'].str[5:7].str.contains('GC')])
    other = other[~other['sequence'].str[5:7].str.contains('GC')]
    print 'AT:'
    print len(df[df['looks like'] == 'AT'])
    print len(other[other['sequence'].str[5:7].str.contains('AT')])
    other = other[~other['sequence'].str[5:7].str.contains('AT')]
    print 'AG:'
    print len(df[df['looks like'] == 'AG'])
    print len(other[other['sequence'].str[3:5].str.contains('AG')])
    other = other[~other['sequence'].str[3:5].str.contains('AG')]
    print 'AC:'
    print len(df[df['looks like'] == 'AC'])
    print len(other[other['sequence'].str[3:5].str.contains('AC')])
    other = other[~other['sequence'].str[3:5].str.contains('AC')]

    print (len(other[other['sequence'].str[6:8].str.contains('AG')]))
    print (len(other[other['sequence'].str[6:8].str.contains('AC')]))
    print (len(other[other['sequence'].str[6:8].str.contains('TT')]))
    print (len(other[other['sequence'].str[6:8].str.contains('TA')]))



    print 'Other sequences:'
    print len(other[other['looks like'] == ''])
    return other
    
def compare_peak_junc_df(peak_df, junc_df, organism = None):
    print str(len(peak_df))+' peaks'
    print str(len(junc_df))+' junctions'
    
    new_df = pd.DataFrame(columns=peak_df.columns)

    junc_type = []
    ann_seq1 = []
    ann_seq2 = []
    junc_size = []
    ann_size = []
    junc_coords = []
    ann_coords = []
    match_count = 0
    for tx in list(set(peak_df['transcript'].tolist())):
        tx_peak = peak_df[peak_df['transcript'] == tx]
        if organism == 'pombe':
            tx_junc = junc_df[junc_df['transcript'] == tx+'.1']
        else:
            tx_junc = junc_df[junc_df['transcript'] == tx+'T0']
        for index, row in tx_peak.iterrows():
            match_flag = False
            #print tx_junc['start']
            peak_range = range(row['position']-1,row['position']+1)
            for pos in peak_range:
                if match_flag is False:
                    tx_juncA = tx_junc[tx_junc['type'] != '5p tethered']
                    tx_juncB = tx_junc[tx_junc['type'] != '3p tethered']
                    
                    if pos in tx_juncA['start'].tolist():
                        match_count += 1
                        match_flag = True
                        for junc_index, junc_row in tx_juncA[tx_juncA['start'] == pos].iterrows():
                            junc_type.append(junc_row['type'])
                            ann_seq1.append(junc_row['annotated sequence1'])
                            ann_seq2.append(junc_row['annotated sequence2'])
                            junc_size.append(junc_row['size'])
                            ann_size.append(junc_row['annotated intron size'])
                            junc_coords.append((junc_row['start'],junc_row['end']))
                            ann_coords.append((junc_row['annotated intron start'],junc_row['annotated intron end']))
                            new_df = new_df.append(row)
                        break
        
                    elif pos in tx_juncB['end'].tolist():
                        match_count += 1
                        match_flag = True
                        for junc_index, junc_row in tx_juncB[tx_juncB['end'] == pos].iterrows():
                            junc_type.append(junc_row['type'])
                            ann_seq1.append(junc_row['annotated sequence1'])
                            ann_seq2.append(junc_row['annotated sequence2'])
                            junc_size.append(junc_row['size'])
                            ann_size.append(junc_row['annotated intron size'])
                            junc_coords.append((junc_row['start'],junc_row['end']))
                            ann_coords.append((junc_row['annotated intron start'],junc_row['annotated intron end']))
                            new_df = new_df.append(row)
                        break                   
    print "Overlap:"
    print match_count
    
    new_df['junction type'] = junc_type
    new_df['annotated sequence1'] = ann_seq1
    new_df['annotated sequence2'] = ann_seq2
    new_df['junction size'] = junc_size
    new_df['annotated intron size'] = ann_size
    new_df['junction coords'] = junc_coords
    new_df['annotated intron coords'] = ann_coords
    return new_df

def peak_seq_enrichment(df, AorT, GorC):
    unpeaks = df[df['type'] == 'other']
    unpeaks = unpeaks.append(df[df['type'] == 'intronic'])
    print "Number of unpredicted peaks:"
    print len(unpeaks)
    nucs = ['G','A','C','T']
    dinucs = set()
    for nuc in nucs:
        for nuc2 in nucs:
            dinucs.add(nuc+nuc2)
    
    five = {}
    three = {}
    for dinuc in dinucs:
        five[dinuc] = len(unpeaks[unpeaks['sequence'].str[6:8].str.contains(dinuc)])
        three[dinuc] = len(unpeaks[unpeaks['sequence'].str[4:6].str.contains(dinuc)])

    p_dict = {'A':AorT, 'T':AorT, 'C':GorC, 'G':GorC}
    five_LO = {}
    three_LO = {}
    for dinuc in five.keys():
        p_dinuc = p_dict[dinuc[0]]*p_dict[dinuc[1]]
        phat_dinuc = five[dinuc]/float(len(unpeaks))
        phat_dinuc2 = three[dinuc]/float(len(unpeaks))

        SE = numpy.sqrt(phat_dinuc*(1-phat_dinuc)/len(unpeaks))
        SE2 = numpy.sqrt(phat_dinuc2*(1-phat_dinuc2)/len(unpeaks))
        Z = (phat_dinuc-p_dinuc)/SE
        Z2 = (phat_dinuc2-p_dinuc)/SE2

        pvalue = stats.norm.sf(Z)
        pvalue2 = stats.norm.sf(Z2)
        LO = numpy.log((1-pvalue)/pvalue)
        LO2 = numpy.log((1-pvalue2)/pvalue2)

        five_LO[dinuc] = LO
        three_LO[dinuc] = LO2

    fig, ax = plt.subplots(figsize=(12,6))
    width = 0.35
    ind = numpy.arange(len(five_LO.keys()))
    rects2 = ax.bar(ind, three_LO.values(), width, color='turquoise', label='Before peak')
    rects1 = ax.bar(ind + width, five_LO.values(), width, color='coral', label='After peak')
    ax.plot([-1,17],[0,0],'-', color='black')
    ax.plot([-1,17],[2.94,2.94], '--', color='0.7', label='95% CI')
    ax.plot([-1,17],[-2.94,-2.94], '--', color='0.7')

    ax.set_xlim([-1,17])
    ax.set_xticklabels(five_LO.keys())
    ax.set_xticks(ind + width / 2)
    ax.set_ylabel('Log odds dinucleotide enrichment')
    ax.set_title('Unpredicted peaks')
    ax.legend()
    
    return fig

def add_intron_size(peaks_df, gff3, organism=None):
    ss_dict, flag = list_splice_sites(gff3, organism=organism)
    ss_dict = SPJunctions.collapse_ss_dict(ss_dict)
    no_peaks = ss_dict
    intron_sizes = []
    for index, row in peaks_df.iterrows():
        if row['type'] != 'intronic':
            intron_sizes.append(numpy.NaN)
        else:
            sites = ss_dict[row['transcript']]
            assigned=False
            for pair in sites:
                if pair[0] > pair[1]:
                    if row['position'] >= pair[1] and row['position'] <= pair[0]:
                        intron_sizes.append(pair[0]-pair[1])
                        assigned=True
                        no_peaks[row['transcript']].remove(pair)
                        break
                else:
                    if row['position'] >= pair[0] and row['position'] <= pair[1]:
                        intron_sizes.append(pair[1]-pair[0])
                        assigned=True
                        no_peaks[row['transcript']].remove(pair)
                        break
            if assigned is False:
                intron_sizes.append(numpy.NaN)
    peaks_df['intron size'] = intron_sizes
    return peaks_df,  no_peaks