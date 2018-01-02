import sys
import os
from subprocess import check_output
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneUtility
import SPTools as SP
from collections import OrderedDict
import csv
import pysam
from matplotlib import pyplot as plt

def count_reads_in_window(bam, chrom, start, end, strand):
    if type(bam) == str:
        bam = pysam.Samfile(bam)
    bam_iter = bam.fetch(chrom, start, end)
    read_count = 0
    for read in bam_iter:
        if strand == "+":
            if read.is_reverse:
                read_count += 1
        elif strand == "-":
            if not read.is_reverse:
                read_count +=1
    return read_count

def tx_info(tx, tx_dict):
    chrom = tx_dict[tx][3]
    strand = tx_dict[tx][2]
    start = tx_dict[tx][0]
    end = tx_dict[tx][1]
    exons = None
    if len(tx_dict[tx]) > 4 and len(tx_dict[tx][4]) > 0:
        CDS_start = min(tx_dict[tx][4])
        CDS_end = max(tx_dict[tx][5])
        exons = zip(tx_dict[tx][4],tx_dict[tx][5])
    else:
        CDS_start = start
        CDS_end = end
        
    return start, end, chrom, strand, CDS_start, CDS_end, exons

def count_aligned_reads(bam_file):
    total = check_output(['samtools','view','-F 0x04','-c',bam_file]).strip()
    total = float(total)/1000000.
    return total

### Builds a pandas series from a bam file for a set of coordinates. Index is genome coordinate and value is number of reads
def generate_read_series(bam_iterator, chrom, start, end, strand, baseline=0):
    s = pd.Series(baseline, index=range(start-25, end+25))
    for read in bam_iterator:
        if read.is_reverse and strand == '+':
            pos = read.reference_start
            if pos not in s.index:
                s[pos] = 0
            s[pos] += 1
        elif not read.is_reverse and strand == '-':
            pos = read.reference_start
            if pos not in s.index:
                s[pos] = 0
            s[pos] = s[pos]+1
    s = s.dropna()
    s = s[s > 0]
    s = s.sort_index()
    return s

### Build a dictionary of read series based on a bam file and transcript dictionary
def map_all_transcripts(gff3, bam_file):
    organism=None
    if 'pombe' in gff3:
        organism = 'pombe'
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    
    bam = pysam.Samfile(bam_file)
    series_dict = {}
    for tx in tx_dict:
        start, end, chrom, strand, CDS_start, CDS_end, exons = tx_info(tx, tx_dict)
        series_dict[tx] = generate_read_series(bam, chrom, start, end, strand)
    return series_dict

def plot_transcripts(series_dict, tx_list):
    if type(tx_list) != list:
        tx_list = [tx_list]
    
    for tx in tx_list:
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(111)
        ax.bar(series_dict[tx].index, series_dict[tx], width, color='darkslateblue')

        plt.show()
        plt.clf()
    
def count_PE_reads(open_bam, chrom, start, end, strand, both_strands=False, count_junctions=False):
    if type(open_bam) != pysam.libcsamfile.Samfile:
        open_bam = pysam.Samfile(open_bam)
    
    iterator = open_bam.fetch(chrom, start-150, end+150)
    
    count = 0
    for read in iterator:
        if both_strands is False:
            # For reads that start or end in the intron
            if read.reference_start in range(start, end) or read.reference_end in range(start, end):
                if not read.is_reverse and strand == '+': count += 1
                elif read.is_reverse and strand == '-': count += 1

            # For reads that span the intron
            elif read.reference_start <= start and read.reference_end >= end:
                intron = False
                if count_junctions is False:
                    # Check to see if the read contains a junction
                    if len(read.cigartuples) == 3 and read.cigartuples[1][0] == 3:
                        intron = True

                if intron is False:
                    if not read.is_reverse and strand == '+': count += 1
                    elif read.is_reverse and strand == '-': count += 1
                            
        else:
            # Don't need to worry about strand or read1 vs. read2, otherwise same as above
            if read.reference_start in range(start, end) or read.reference_end in range(start, end):
                count += 1
            elif read.reference_start <= start and read.reference_end >= end:
                intron = False
                if count_junctions is False:
                    # Check to see if the read contains a junction
                    if len(read.cigartuples) == 3 and read.cigartuples[1][0] == 3:
                        intron = True
                if intron is False:
                    count += 1
                    
    return count

def PE_intron_retention_from_annotation(bam_list, organism, both_strands=False, count_junctions=False):
    if 'crypto' in organism.lower():
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
        organism=None
    elif 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        organism='pombe'
    elif 'cerev' in organism.lower():
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        organism=None
        
    tx_dict = GT.build_transcript_dict(gff3, organism=organism)
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss_dict)
    #ss_dict = {k:v for k, v in ss_dict.items() if k in ss_dict.keys()[:100]}
    
    open_bams = {}
    data_dict = {}
    for bam in bam_list:
        open_bams[bam] = pysam.Samfile(bam)
        name = bam.split('/')[-1].split('_sorted.bam')[0]
        print name
        data_dict[bam] = {bam:name, 'reads in transcript':[], 'reads in intron':[]}
        
    column_dict = {'transcript':[],'intron start':[],'intron end':[],'transcript size':[],'chromosome':[],'strand':[]}
    
    for tx, splice_sites in ss_dict.iteritems():
        print tx
        # Get information for overall transcript
        if organism == 'pombe': iso = tx+'.1'
        else: iso = tx+'T0'
        start, end, chrom, strand, CDS_start, CDS_end, exons = GT.tx_info(iso, tx_dict)

        tx_counts = {}
        for bam, open_bam in open_bams.iteritems():
            # Count reads in transcript
            tx_counts[bam] = count_PE_reads(open_bam, chrom, start, end, strand, both_strands=both_strands, count_junctions=True)
        
        # Iterate over all annotated introns in transcript
        for five, three in splice_sites:
            column_dict['transcript'].append(tx)
            column_dict['intron start'].append(five)
            column_dict['intron end'].append(three)
            column_dict['transcript size'].append(end-start)
            column_dict['chromosome'].append(chrom)
            column_dict['strand'].append(strand)
            
            for bam, open_bam in open_bams.iteritems():
                if strand == '+':
                    intron_counts = count_PE_reads(open_bam, chrom, five, three, strand, both_strands=both_strands, count_junctions=count_junctions)
                if strand == '-':
                    intron_counts = count_PE_reads(open_bam, chrom, three, five, strand, both_strands=both_strands, count_junctions=count_junctions)
                
                data_dict[bam]['reads in transcript'].append(tx_counts[bam])
                data_dict[bam]['reads in intron'].append(intron_counts)
                
    df = pd.DataFrame(columns=column_dict.keys(), index=range(len(column_dict['transcript'])))
    for col, info in column_dict.iteritems():
        df[col] = info
    df['intron size'] = (df['intron start']-df['intron end']).apply(abs)
    
    for bam, data in data_dict.iteritems():
        df[data[bam]+': reads in transcript'] = data['reads in transcript']
        df[data[bam]+': reads in intron'] = data['reads in intron']
        df[data[bam]+': intron retention'] = ((df[data[bam]+': reads in intron']/float(sum(df[data[bam]+': reads in intron'])))/
                                              (df[data[bam]+': reads in transcript']/float(sum(df[data[bam]+': reads in transcript']))))
        
    return df