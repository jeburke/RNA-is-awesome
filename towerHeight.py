__author__ = 'jordanburke'

'''Usage: python towerHeight.py <exon_table> <total_table> <transcipt_lengths> <configuration file> <RNAi_list> <prefix>

Configuration file format (tab separated):
SampleName1 SampleName2
ControlRNAcountsA   ControlRNAcountsA
ControlRNAcountsB   ControlRNAcountsB

Note: .bam files must be named *A_al_sorted.bam for exon sequencing and *B_al_sorted.bam for total sequencing'''

import pandas
import numpy
import matplotlib.pyplot as plt
import sys
from scipy import stats
import math
import beeswarm
import random
import argparse

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
## Build dataframes from tables written by analyzeSp.py. Reads ##
## file containing Transcript column and columns with counts   ##
## for each type of read (exon, intron or total)               ##
#################################################################

def build_tables(file):
    fin = open(file, "r")
    df = pandas.read_table(fin)
    df.columns = df.columns.str.rstrip("_al.sorted.bam")
    df.fillna(0)
    return df

#################################################################
## Normalize counts in A samples to counts in B samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA.                                       ##
#################################################################

#def normalize_counts(table1, table2, table3, c1, c2, c3, c4):
def normalize_counts(table1,table2,table3,config):
    df1 = table1
    df2 = table2
    df3 = table3
    df4 = config
    merged = pandas.merge(df1,df2,on="Transcrip",left_index=True,how='left',suffixes=('_exons','_total'))
    merged = pandas.merge(merged,df3, on="Transcrip", left_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df4.columns:
        configNames.append(name)
        s = pandas.Series(df4[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    for name in merged.columns:
        names.append(name)
        if name.strip() == "Transcrip":
            print "Reading table"
        elif n < len(configNames) and name == configNames[n]+"-A_exons":
            c1 = controlValues[n][0]
            print c1
            c2 = controlValues[n][1]
            print c2
            merged[name+"_norm"] = pandas.Series((merged[name]/c1)/(merged[name.strip("A_exons")+"B_total"]/merged[name.strip("Length")]/c2), index = merged.index)
            n += 1
    merged.fillna(0)
    return merged

##################################################################
## Filter transcripts based on an input table. In this case I'm ##
## giving it a list of RNAi targets. Can be anything, 1 column  ##
## format with CNAG IDs                                         ##
##################################################################

def filter_transcripts(mergedtable, list):
    geneID = []
    fin = open(list, "r")
    df = pandas.DataFrame(mergedtable)
    df = df.set_index('Transcrip')
    for line in fin:
        if line.startswith("CNAG"):
            gene = line.strip()
            geneID.append(gene)
            geneID.append(gene+"T0")
            geneID.append(gene+"T1")
            geneID.append(gene+"T2")
            geneID.append(gene+"T3")
            geneID.append(gene+"T4")
            geneID.append(gene+"T5")
    RNAi_df = df[(df.index.isin(geneID))]
    return RNAi_df

#################################################################
## Convert output from normalize_counts to lists for plotting  ##
#################################################################

def get_ratios(normalizedResults, a):
    name = ""
    values = []
    n = 0
    for name in normalizedResults.columns:
        if name[-5:] == "_norm":
            if n == a:
                s = pandas.Series(normalizedResults[name])
                values = s.tolist()
                print name
                n += 1
            else:
                n += 1
    return values


#################################################################
## Convert input tables to dataframes                          ##
#################################################################

exonCounts = build_tables(sys.argv[1])
totalCounts = build_tables(sys.argv[2])
transcriptLength = build_tables(sys.argv[3])
configFile = build_tables(sys.argv[4])

#################################################################
## Process exon counts and filter for genes of interest        ##
#################################################################

#normalizedTable = normalize_counts(exonCounts,totalCounts, transcriptLength, int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]))
normalizedTable = normalize_counts(exonCounts,totalCounts, transcriptLength, configFile)
filteredList = filter_transcripts(normalizedTable, sys.argv[5])

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

xvalues = get_ratios(normalizedTable, 0)
xvalues = [0 if math.isnan(x) else x for x in xvalues]
yvalues = get_ratios(normalizedTable, 1)
yvalues = [0 if math.isnan(x) else x for x in yvalues]

RNAixvalues = get_ratios(filteredList, 0)
RNAixvalues = [0 if math.isnan(x) else x for x in RNAixvalues]
RNAiyvalues = get_ratios(filteredList, 1)
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
ax1.set_xlim([-1,90])
ax1.set_ylim([-1,90])
plt.show()

#randm_xvalues = [xvalues[i] for i in sorted(random.sample(xrange(len(xvalues)), 100))]
#randm_RNAixvalues = [RNAixvalues[i] for i in sorted(random.sample(xrange(len(RNAixvalues)), 100))]
#randm_yvalues = [yvalues[i] for i in sorted(random.sample(xrange(len(yvalues)), 100))]
#randm_RNAiyvalues = [RNAiyvalues[i] for i in sorted(random.sample(xrange(len(RNAiyvalues)), 100))]
#bs, ax2 = beeswarm.beeswarm([randm_xvalues,randm_RNAixvalues, randm_yvalues,randm_RNAiyvalues], positions = [1,2,3,4], labels = ["Replicate 1 All", "Replicate 1 RNAi","Replicate 2 All","Replicate 2 RNAi"],col=["blue","coral","blue","coral"])
#plt.show()
