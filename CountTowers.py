__author__ = 'jordanburke'
import HTseq
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
import gobs
import numpy
import pandas as pd
import copy
import os
import multiprocessing
import itertools
import time
import shelve
import collections
import analyzeSplicing

BAM_files=[]
BAM_file_readers=[]
coverage_vectors=[]
directories=[]
standard_values=[]

#Read configuration file with all file locations

config_file=open(sys.argv[1], "r")
for line in config_file:
    data=line.split("=")
    if data[0]=="GFF":
        GFF_file_name=data[1].strip()
        print GFF_file_name
    if data[0]=="BAM":
        BAM_file_name=data[1].strip()
        print BAM_file_name
        BAM_files.append(BAM_file_name)
        directoryname=basedir+BAM_file_name.split("/")[-1]+"results"
        if os.path.exists(directoryname)==False:
            os.mkdir(directoryname)
            os.mkdir(directoryname+"/wigs")
            os.mkdir(directoryname+"/spreadsheets")
            os.mkdir(directoryname+"/plots")
        directories.append(directoryname)
    if data[0]=="STANDARD":
        standard_values=data[1].split(",")
        print standard_values

config_file.close()
print 'done'

#Build genome database

print "building H99 genome database..."
genome,chromosomes=gobs.buildGenomeDatabase(basedir+"conf/H99.config")

t2=time.time()

print "Execution time=", int(t2-t1), " seconds"

#Count reads that end at each 5' splice site and normalize to standard reads

transcript_feature_di={}
CDS_feature_di={}
transcript_len_di={}
CDS_len_di={}
exon_feature_di={}

for chromosome in chromosomes:
    for gene in chromosome.getElements():
        for transcript in gene.getTranscripts():
            transcript_coord=transcript.getCoordinates()
            CDS_coord=transcript.getCDSCoordinates()
            transcript_iv=gobs.HTSeq_iv(transcript_coord)
            CDS_iv=gobs.HTSeq_iv(CDS_coord)
            transcript_feature_di[transcript.getName()]=transcript_iv
            transcript_len_di[transcript.getName()]=transcript_iv.end-transcript_iv.start
            CDS_feature_di[transcript.getName()]=CDS_iv
            CDS_len_di[transcript.getName()]=CDS_iv.end-CDS_iv.start

#            exon_coord=exon_starts_stops.getCoordinates()
#            exon_iv=gobs.HTSeq_iv(exon_coord)
#            exon_feature_di[transcript.getName()]=exon_iv

print len(transcript_feature_di), " transcript intervals defined"
print len(CDS_feature_di)," CDS feature intervals defined"

df=pd.DataFrame()
column_headers=[]

p=multiprocessing.Pool(processes=4)

params=[]

for name in BAM_files:
    column_headers.append(name.split("/")[-1])

    if "B" in name:
        params.append((name,transcript_feature_di))

BResults=p.map(gobs.count_reads,params)

#for name in BAM_files:
#   column_heads.append(name.split("/")[-1])

#    if "A" in name:
#       params.append((name,exon_feature_di))

#AResults=p.map(analyzeSplicing.count_reads_end,params)

for result, column_header in zip(results,column_headers):
    ds=pd.Series(result)
    df[column_header]=ds

transcript_len_series=pd.Series(transcript_len_di)
CDS_len_series=pd.Series(CDS_len_di)
df["transcript_length"]=transcript_len_series
df["CDS_length"]=CDS_len_series



df.to_csv(basedir+"/counts.txt",sep="\t")