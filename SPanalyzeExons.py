__author__ = 'jordanburke'

import HTSeq
import matplotlib
matplotlib.use('TkAgg')
import gobsChr
import pandas as pd
import copy
import os
import multiprocessing
import itertools
import time
import shelve
import collections
import sys

t1=time.time()
print 'reading configuration file'
BAM_files=[]
BAM_file_readers=[]
coverage_vectors=[]
directories=[]

#############################################
#  This code reads a configuration file     #
#  that contains the paths to the .bam files#
#  to be analyzed and the GFF3 file used to #
#  generate those bamfiles using TOPHAT     #
#############################################

config_file=open(sys.argv[1])
for line in config_file:
    data=line.split("=")
    if data[0]=="GFF3":
        GFF_file_name=data[1].strip()
        print GFF_file_name
    if data[0]=="BAM":
        BAM_file_name=data[1].strip()
        print BAM_file_name
        BAM_files.append(BAM_file_name)
config_file.close()
for f in BAM_files:
    BAM_file_readers.append(HTSeq.BAM_Reader(f))

######################################################################################
#                                                                                    #
#   One line of code to build a transcript object generator                          #
#                                                                                    #
######################################################################################

transcripts=gobsChr.transcript_factory(sys.argv[1])

###################################################################################################################
# COUNTING READS                                                                                                  #
# First name a dictionaries of features of interest in which the keys are names and the values are HT-Seq Genomic #
# Interval objects.  Then use gobs.count_read which returns a collections.Counter objects of names and counts     #
###################################################################################################################

transcript_feature_di={}
transcript_len_di={}
CDS_feature_di={}
CDS_len_di={}

for transcript in transcripts:
    transcript_coord=transcript.getCoordinates()
    transcript_iv=gobsChr.HTSeq_iv(transcript_coord)
    transcript_feature_di[transcript.getName()]=transcript_iv
    transcript_len_di[transcript.getName()]=transcript_iv.end-transcript_iv.start
    transcript_exon_coords=transcript.getAllExonCoordinates()
    transcript_coord=transcript.getCoordinates()
    CDS_coord=transcript.getCDSCoordinates()
    CDS_iv=gobsChr.HTSeq_iv(CDS_coord)
    CDS_feature_di[transcript.getName()]=CDS_iv
    CDS_len_di[transcript.getName()]=CDS_iv.end-CDS_iv.start

print len(transcript_feature_di), " transcript intervals defined"
print len(CDS_feature_di)," CDS feature intervals defined"

df_total=pd.DataFrame()
column_headers=[]
p=multiprocessing.Pool(processes=24)
params_total_counts=[]

for bamfilereader in BAM_file_readers:
    column_headers.append(bamfilereader.filename.split("/")[-1])
    params_total_counts.append((bamfilereader,CDS_feature_di))

total_results=p.map(gobsChr.count_reads,params_total_counts)

for result, column_header in zip(total_results,column_headers):
    ds=pd.Series(result)
    df_total[column_header]=ds

df_total_sorted=df_total.sort()

fout = open("{0}_totalcounts.tsv".format(sys.argv[2]), "w")
fout.write(pd.DataFrame.to_csv(df_total_sorted, sep='\t'))
