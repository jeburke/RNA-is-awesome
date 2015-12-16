__author__ = 'jordanburke'

'''This script compares counts for exons in SP-B samples to total reads for that transcript in the mature mRNA samples (SP-W samples)
This was the closest equivalent to the "Spliceosome occupancy" assay from Dumesic et al. 2013.

Usage: python NormalizeToMature.py <SPanalyzeExons output> <transcipt_lengths> <configuration file> <RNAi_list> <prefix>

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
import random

def configure(file):
    sample = []
    controlReads = []
    fin = open(file, "r")
    for line in fin:
        row = line.split("\t")
        sample.append(row[0])
        controlReads.append(row[1])

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

#################################################################
## Get lists from dataframes to make scatter plots and do some ##
## statistics                                                  ##
#################################################################

xvalues = SPTools.get_ratios(normalizedTable, 0)
xvalues = [0 if math.isnan(x) else x for x in xvalues]
yvalues = SPTools.get_ratios(normalizedTable, 1)
yvalues = [0 if math.isnan(x) else x for x in yvalues]

RNAixvalues = SPTools.get_ratios(filteredList, 0)
RNAixvalues = [0 if math.isnan(x) else x for x in RNAixvalues]
RNAiyvalues = SPTools.get_ratios(filteredList, 1)
RNAiyvalues = [0 if math.isnan(x) else x for x in RNAiyvalues]

#################################################################
## Histogram                                                   ##
#################################################################

fig = plt.figure()
ax1 = fig.add_subplot(111)
#ax1.hist(xvalues, bins = 60, normed=True, color='royalblue', label='All')
#ax1.hist(RNAixvalues, bins = 20, normed=True, color='coral', alpha=0.5, label='RNAi')
#ax1.set_xlabel("Value")
#ax1.set_ylabel("Frequency")
ax1.scatter(xvalues, yvalues, alpha = 0.5)
#ax1.set_xscale("log")
#ax1.set_yscale("log")
ax1.set_xlim(0,0.2)
ax1.set_ylim(0,0.2)
#ax1.legend()
plt.show()