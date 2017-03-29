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

def count_aligned_reads(bam_file):
    total = check_output(['samtools','view','-F 0x04','-c',bam_file]).strip()
    total = float(total)/1000000.
    return total