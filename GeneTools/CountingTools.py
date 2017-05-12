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
    
    