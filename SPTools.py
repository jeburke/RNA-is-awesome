__author__ = 'jordanburke'

import pandas
import numpy
import matplotlib.pyplot as plt
import sys
from scipy import stats
import math
import beeswarm
import random
import argparse


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
            merged[name+"_norm"] = pandas.Series((merged[name]/c1)/(merged[name.strip("A_exons")+"B_total"]/c2), index = merged.index)
            n += 1
    merged.fillna(0)
    return merged


#################################################################
## Normalize counts in B samples to counts in W samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA.                                       ##
#################################################################

def normalize_to_mature(table1, table2, table3, config):
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
        elif n < len(configNames) and name == configNames[n]+"-B_total":
            c1 = controlValues[n][0]
            print c1
            c2 = controlValues[n][1]
            print c2
            c3 = controlValues[n][2]
            print c3
            c4 = controlValues[n][3]
            print c4
            merged[name+"_norm"] = pandas.Series((merged[name]/c3)/(merged["150831CM763-W_total"]/c4/merged[name.strip("Length")]), index = merged.index)
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