__author__ = 'jordanburke'

'''This script compares counts for exons in SP-B samples to total reads for that transcript in the mature mRNA samples (SP-W samples)
This was the closest equivalent to the "Spliceosome occupancy" assay from Dumesic et al. 2013.

Usage: python NormalizeToMature.py <cleaved 5p exon counts> <cleaved intron counts> <premRNA total counts> <mature RNA counts> <transcipt_lengths> <configuration file> <RNAi_list> <prefix>

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

#################################################################
## Convert input tables to dataframes                          ##
#################################################################

ExonCounts = SPTools.build_tables(sys.argv[1])
IntronCounts = SPTools.build_tables(sys.argv[2])
SPtotalCounts = SPTools.build_tables(sys.argv[3])
totalCounts = SPTools.build_tables(sys.argv[4])

transcriptLength = SPTools.build_tables(sys.argv[5])
configFile = SPTools.build_tables(sys.argv[6])

#################################################################
## Process exon counts and filter for genes of interest        ##
#################################################################

SPvM = SPTools.normalize_to_mature(SPtotalCounts, totalCounts, transcriptLength, configFile, "B")
SPvM_RNAi = SPTools.filter_transcripts_by_cnag(SPvM, sys.argv[7])
print SPvM

ExvM = SPTools.normalize_to_mature(ExonCounts, totalCounts, transcriptLength, configFile, "A")
ExvM_RNAi = SPTools.filter_transcripts_by_cnag(ExvM, sys.argv[7])

IntvM = SPTools.normalize_to_mature(IntronCounts, totalCounts, transcriptLength, configFile, "A")
IntvM_RNAi = SPTools.filter_transcripts_by_cnag(IntvM, sys.argv[7])

#################################################################
## Write 2 files - one with all data and one with filtered set ##
#################################################################

fout = open("{0}_occupancy.tsv".format(sys.argv[5]), "w")
fout.write(pandas.DataFrame.to_csv(SPvM, sep='\t'))

fout = open("{0}_RNAi_occupancy.tsv".format(sys.argv[5]),"w")
fout.write(pandas.DataFrame.to_csv(SPvM_RNAi, sep='\t'))

#################################################################
## Get lists from dataframes to make scatter plots and do some ##
## statistics                                                  ##
#################################################################

SPvM_xvalues = SPTools.get_ratios(SPvM, 0)
SPvM_yvalues = SPTools.get_ratios(SPvM, 1)
SPvM_RNAixvalues = SPTools.get_ratios(SPvM_RNAi, 0)
SPvM_RNAiyvalues = SPTools.get_ratios(SPvM_RNAi, 1)

ExvM_xvalues = SPTools.get_ratios(ExvM, 0)
ExVM_yvalues = SPTools.get_ratios(ExvM, 1)
ExvM_RNAixvalues = SPTools.get_ratios(ExvM_RNAi, 0)
ExvM_RNAiyvalues = SPTools.get_ratios(ExvM_RNAi, 1)

InvM_xvalues = SPTools.get_ratios(IntvM, 0)
InvM_yvalues = SPTools.get_ratios(IntvM, 1)
InvM_RNAixvalues = SPTools.get_ratios(IntvM_RNAi, 0)
InvM_RNAiyvalues = SPTools.get_ratios(IntvM_RNAi, 1)

logSPvM_xvalues = SPTools.log_ratios(SPvM_xvalues)
logSPvM_yvalues = SPTools.log_ratios(SPvM_yvalues)
logSPvM_RNAixvalues = SPTools.log_ratios(SPvM_RNAixvalues)
logSPvM_RNAiyvalues = SPTools.log_ratios(SPvM_RNAiyvalues)

SPTools.scatter_plot(logSPvM_xvalues,logSPvM_yvalues,logSPvM_RNAixvalues,logSPvM_RNAiyvalues, "SP total/Mature RNA", "All", "RNAi targets")


