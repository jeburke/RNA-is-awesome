__author__ = 'jordanburke'
__author__ = 'jordanburke'

'''This script compares counts at splice site features in SP-A to total reads for that transcript in the total Prp19 IP samples (SP-B samples)

Usage: python NormalizeToMature.py <cleaved 5p exon counts> <cleaved intron counts> <premRNA total counts> <transcipt_lengths> <configuration file> <RNAi_list> <prefix>

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

exonCounts = SPTools.build_tables(sys.argv[1])
intronCounts = SPTools.build_tables(sys.argv[2])
totalCounts = SPTools.build_tables(sys.argv[3])
transcriptLength = SPTools.build_tables(sys.argv[4])
configFile = SPTools.build_tables(sys.argv[5])

#################################################################
## Process exon counts and filter for genes of interest        ##
#################################################################

normalizedExonTable = SPTools.normalize_AtoB(exonCounts, totalCounts, transcriptLength, configFile)
filteredExonList = SPTools.normalize_AtoB_filter_by_counts(exonCounts, totalCounts, transcriptLength, configFile, 1000)
RNAiExons = SPTools.filter_transcripts_by_cnag(normalizedExonTable, sys.argv[6])

normalizedIntronTable = SPTools.normalize_AtoB(intronCounts, totalCounts, transcriptLength, configFile)
filteredIntronList = SPTools.normalize_AtoB_filter_by_counts(intronCounts, totalCounts, transcriptLength, configFile, 1000)
RNAiIntrons = SPTools.filter_transcripts_by_cnag(normalizedIntronTable, sys.argv[6])

#################################################################
## Write 2 files - one with all data and one with filtered set ##
#################################################################

fout = open("{0}_normalized5pExons.tsv".format(sys.argv[7]), "w")
fout.write(pandas.DataFrame.to_csv(normalizedExonTable, sep='\t'))

fout = open("{0}_normalizedIntrons.tsv".format(sys.argv[7]), "w")
fout.write(pandas.DataFrame.to_csv(normalizedIntronTable, sep='\t'))

fout = open("{0}_RNAi_normalized5pExons.tsv".format(sys.argv[7]),"w")
fout.write(pandas.DataFrame.to_csv(RNAiExons, sep='\t'))

fout = open("{0}_RNAi_normalizedIntrons.tsv".format(sys.argv[7]),"w")
fout.write(pandas.DataFrame.to_csv(RNAiIntrons, sep='\t'))

####################################################################
## Get lists from dataframes to make scatter plots and do some    ##
## statistics - second argument is column number (indexed from 0) ##
####################################################################

Exon_xvalues = SPTools.get_ratios(normalizedExonTable, 0)
Exon_yvalues = SPTools.get_ratios(normalizedExonTable, 1)

Filt_exon_xvalues = SPTools.get_ratios(filteredExonList, 0)
Filt_exon_yvalues = SPTools.get_ratios(filteredExonList, 1)

Intron_xvalues = SPTools.get_ratios(normalizedIntronTable, 0)
Intron_yvalues = SPTools.get_ratios(normalizedIntronTable, 1)

Filt_intron_xvalues = SPTools.get_ratios(filteredIntronList, 0)
Filt_intron_yvalues = SPTools.get_ratios(filteredIntronList, 1)

RNAi_ex_xvalues = SPTools.get_ratios(RNAiExons, 0)
RNAi_ex_yvalues = SPTools.get_ratios(RNAiExons, 1)
RNAi_int_xvalues = SPTools.get_ratios(RNAiIntrons, 0)
RNAi_int_yvalues = SPTools.get_ratios(RNAiIntrons, 1)


#################################################################
## Scatter plot of replicates                                  ##
#################################################################

logExon_xvalue = SPTools.log_ratios(Exon_xvalues)
logExon_yvalue = SPTools.log_ratios(Exon_yvalues)
logFiltex_xvalue = SPTools.log_ratios(Filt_exon_xvalues)
logFiltex_yvalue = SPTools.log_ratios(Filt_exon_yvalues)
logRNAi_ex_xvalue = SPTools.log_ratios(RNAi_ex_xvalues)
logRNAi_ex_yvalue = SPTools.log_ratios(RNAi_ex_yvalues)
print "Number of 5' exons: "+str(len(logExon_xvalue))
print "Number of 5' exons above cutoff: "+str(len(logFiltex_xvalue))
SPTools.scatter_plot(logExon_xvalue,logExon_yvalue,logFiltex_xvalue,logFiltex_yvalue)
SPTools.scatter_plot(logExon_xvalue,logExon_yvalue,logRNAi_ex_xvalue,logRNAi_ex_yvalue)

logIntron_xvalue = SPTools.log_ratios(Intron_xvalues)
logIntron_yvalue = SPTools.log_ratios(Intron_yvalues)
logFiltin_xvalue = SPTools.log_ratios(Filt_intron_xvalues)
logFiltin_yvalue = SPTools.log_ratios(Filt_intron_yvalues)
print "Number of introns: "+str(len(logIntron_xvalue))
print "Number of introns above cutoff: "+str(len(logFiltin_xvalue))
SPTools.scatter_plot(logIntron_xvalue,logIntron_yvalue,logFiltin_xvalue,logFiltin_yvalue)