__author__ = 'jordanburke'
import HTSeq
import matplotlib
matplotlib.use('TkAgg')
import gobsChr
import pandas as pd
import multiprocessing
import time
import sys

if len(sys.argv) < 2:
    print '''Usage: python analyzeSP.py <configuration file> <output prefix>
             Configuration file format:
             GFF3=/Users/jordanburke/GENOMES/CNA3_FINAL_CALLGENES_1_Chr.gff3
             FASTA=/Users/jordanburke/GENOMES/CNA3-Chr.fa
             CEN=/Users/jordanburke/python-hm/H99centromeres.txt
             SPECIES=CRYPTOCOCCUSNEOFORMANSVARGRUBII
             BAM=/Users/jordanburke/Documents/JEB/HTseq/JEB003/BAMFILES/CM763-A_sorted.bam'''
    sys.exit()

t1=time.time()
#basedir="./"
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
    if data[0]=="GFF":
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

print transcripts

###################################################################################################################
# COUNTING READS                                                                                                  #
# First name a dictionary of features of interest in which the keys are names and the values are HT-Seq Genomic #
# Interval objects.  Then use gobs.count_read which returns a collections.Counter objects of names and counts     #
###################################################################################################################
transcript_feature_di={}
transcript_len_di={}
exon_di={}


for transcript, genome in transcripts:
    transcript_coord=transcript.getCoordinates()
    transcript_iv=gobsChr.HTSeq_iv(transcript_coord)
    transcript_feature_di[transcript.getName()]=transcript_iv
    transcript_len_di[transcript.getName()]=transcript_iv.end-transcript_iv.start
    transcript_exon_coords=transcript.getAllExonCoordinates()
    if len(transcript_exon_coords)==0:
        print "no introns in ",transcript
        continue
    if transcript_exon_coords[0].getStrand() == "Watson":
        transcript_exon_coords.sort()
    else:
        transcript_exon_coords.sort(reverse=True)
    transcript_exon_ivs=[gobsChr.HTSeq_iv(x) for x in transcript_exon_coords]
    #print transcript_exon_ivs
    flag=False
    for interval in transcript_exon_ivs:
        if interval.end-interval.start<1:
            print "transcript", transcript," has a problem with one or more of its exons:"
            print transcript_exon_coords
            flag=True
    if flag==False: exon_di[transcript.getName()]=transcript_exon_ivs

t2=time.time()
print "Execution time=", int(t2-t1), " seconds"

print len(transcript_feature_di), " transcript intervals defined"
print len(exon_di)," exon feature intervals defined"

df_mature = pd.DataFrame()

column_headers=[]
p=multiprocessing.Pool(processes=2)


params_cleaved_exon_counts=[]
df_cleaved_exon_counts=pd.DataFrame()

params_mature_counts = []
for bamfilereader in BAM_file_readers:
    params_mature_counts.append((bamfilereader,exon_di))
    params_cleaved_exon_counts.append((bamfilereader,exon_di))

cleaved_exon_count_results=p.map(gobsChr.count_splicing_intermediates,params_cleaved_exon_counts)
#print cleaved_exon_count_results


mature_count_results=p.map(gobsChr.count_reads_per_region,params_mature_counts)
#print mature_count_results


print len(cleaved_exon_count_results)
print len(mature_count_results), "mRNA Count Results Produced"

for result, column_header in zip(mature_count_results,column_headers):
    df_mature.index.name = "Transcript"
    ds=pd.Series(result)
    df_mature[column_header]=ds

for result,column_header in zip(cleaved_exon_count_results,column_headers):
    index=pd.MultiIndex.from_tuples(result.keys(),names=["Transcript","Intron"])
    ds=pd.Series(result,index=index)
    df_cleaved_exon_counts[column_header]=ds

df_mature_sorted=df_mature.sort()
print df_mature_sorted
df_cleaved_exon_counts_sorted = df_cleaved_exon_counts.sort()

df_mature_sorted.to_csv("{0}_mRNAcounts.tsv".format(sys.argv[2]),sep="\t")
df_cleaved_exon_counts_sorted.to_csv("{0}_exoncounts.tsv".format(sys.argv[2]),sep="\t")