__author__ = 'jordanburke'

'''Usage: python NormalizeToMature.py <SPanalyzeExons output> <transcipt_lengths> <configuration file> <RNAi_list> <prefix>

Configuration file format (tab separated):
SampleName1 SampleName2
ControlRNAcountsA   ControlRNAcountsA
ControlRNAcountsB   ControlRNAcountsB
ControlRNAcountsW   ControlRNAcountsW

Note: .bam files must be named *A_al_sorted.bam for exon sequencing, *B_al_sorted.bam for total spliceosome sequencing, *W_al_sorted.bam for mature sequencing'''

import pandas
import matplotlib.pyplot as plt
import sys
from scipy import stats
import math
import SPTools

def configure(file):
    sample = []
    controlReads = []
    fin = open(file, "r")
    for line in fin:
        row = line.split("\t")
        sample.append(row[0])
        controlReads.append(row[1])
    print sample
    print controlReads

#################################################################
## Convert input tables to dataframes                          ##
#################################################################

totalCounts = SPTools.build_tables(sys.argv[1])
transcriptLength = SPTools.build_tables(sys.argv[2])
configFile = SPTools.build_tables(sys.argv[3])

#################################################################
## Process exon counts and filter for genes of interest        ##
#################################################################

normalizedTable = SPTools.normalize_to_mature(totalCounts, transcriptLength, configFile)
filteredList = SPTools.filter_transcripts(normalizedTable, sys.argv[4])

#################################################################
## Write 2 files - one with all data and one with filtered set ##
#################################################################

fout = open("{0}_occupancy.tsv".format(sys.argv[5]), "w")
fout.write(pandas.DataFrame.to_csv(normalizedTable, sep='\t'))

fout = open("{0}_RNAi_occupancy.tsv".format(sys.argv[5]),"w")
fout.write(pandas.DataFrame.to_csv(filteredList, sep='\t'))
