__author__ = 'jordanburke'

'''Usage: python NormalizeToMature.py <exon_table> <total_table> <transcipt_lengths> <configuration file> <RNAi_list> <prefix>

Configuration file format (tab separated):
SampleName1 SampleName2
ControlRNAcountsA   ControlRNAcountsA
ControlRNAcountsB   ControlRNAcountsB
TotalcountsW   TotalcountsW

Note: .bam files must be named *A_al_sorted.bam for exon sequencing, *B_al_sorted.bam for total spliceosome sequencing, *W_al_sorted.bam for mature sequencing'''

import pandas
import numpy
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

exonCounts = SPTools.build_tables(sys.argv[1])
totalCounts = SPTools.build_tables(sys.argv[2])
transcriptLength = SPTools.build_tables(sys.argv[3])
configFile = SPTools.build_tables(sys.argv[4])

#################################################################
## Process exon counts and filter for genes of interest        ##
#################################################################

normalizedTable = SPTools.normalize_to_mature(exonCounts,totalCounts, transcriptLength, configFile)
filteredList = SPTools.filter_transcripts(normalizedTable, sys.argv[5])

#################################################################
## Write 2 files - one with all data and one with filtered set ##
#################################################################

fout = open("{0}_towers.tsv".format(sys.argv[6]), "w")
fout.write(pandas.DataFrame.to_csv(normalizedTable, sep='\t'))

fout = open("{0}_RNAi.tsv".format(sys.argv[6]),"w")
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
## Two tailed T-test                                           ##
#################################################################

replicate1ttest = stats.ttest_ind(xvalues,RNAixvalues)
print "T-statistic and p-value for replicate 1"
print replicate1ttest
replicate2ttest = stats.ttest_ind(yvalues,RNAiyvalues)
print "T-statistic and p-value for replicate 2"
print replicate2ttest

replicate1mwtest = stats.mannwhitneyu(xvalues,RNAixvalues)
print "Mann-Whitney U test p-value for replicate 1"
print str(replicate1mwtest)
replicate2mwtest = stats.mannwhitneyu(yvalues,RNAiyvalues)
print "Mann-Whitney U test p-value for replicate 2"
print str(replicate2mwtest)

#################################################################
## Scatter plot                                                ##
#################################################################

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(xvalues, yvalues, color = 'royalblue', label="All")
ax1.scatter(RNAixvalues, RNAiyvalues, color = 'coral', label="RNAi targets")
ax1.set_xlabel("Replicate 1 Exons/Total")
ax1.set_ylabel("Replicate 2 Exons/Total")
ax1.legend(loc = 'lower right')
#ax1.set_xlim([-1,90])
#ax1.set_ylim([-1,90])
plt.show()