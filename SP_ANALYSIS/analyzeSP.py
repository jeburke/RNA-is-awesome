import HTSeq
import matplotlib
matplotlib.use('TkAgg')
import gobsChr
import pandas as pd
import multiprocessing
import time
import sys

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

###################################################################################################################
# COUNTING READS                                                                                                  #
# First name a dictionary of features of interest in which the keys are names and the values are HT-Seq Genomic #
# Interval objects.  Then use gobs.count_read which returns a collections.Counter objects of names and counts     #
###################################################################################################################

transcript_feature_di={}
transcript_len_di={}
exons_di={}
introns_di={}
CDS_feature_di={}
CDS_len_di={}

for transcript, genome in transcripts:
    print transcript
    transcript_coord=transcript.getCoordinates()
    transcript_iv=gobsChr.HTSeq_iv(transcript_coord)
    transcript_feature_di[transcript.getName()]=transcript_iv
    transcript_len_di[transcript.getName()]=transcript_iv.end-transcript_iv.start
    transcript_exon_coords=transcript.getAllExonCoordinates()
    transcript_intron_coords=transcript.getAllIntronCoordinates()
    CDS_coord=transcript.getCDSCoordinates()
    CDS_iv=gobsChr.HTSeq_iv(CDS_coord)
    CDS_feature_di[transcript.getName()]=CDS_iv
    CDS_len_di[transcript.getName()]=CDS_iv.end-CDS_iv.start
    if len(transcript_exon_coords)==0:
        print "no introns in ",transcript
        continue
    if transcript_exon_coords[0].getStrand() == "Watson":
        transcript_exon_coords.sort()
    else:
        transcript_exon_coords.sort(reverse=True)
    transcript_exon_ivs=[gobsChr.HTSeq_iv(x) for x in transcript_exon_coords]
    flag=False
    for interval in transcript_exon_ivs:
        if interval.end-interval.start<1:
            print "transcript", transcript," has a problem with one or more of its exons:"
            print transcript_exon_coords
            flag=True
    if flag==False: exons_di[transcript.getName()]=transcript_exon_ivs

    if len(transcript_intron_coords)==0:
        print "no introns in ", transcript
        continue
    if transcript_intron_coords[0].getStrand() == "Watson":
        transcript_intron_coords.sort()
    else:
        transcript_intron_coords.sort(reverse=True)
    transcript_intron_ivs = [gobsChr.HTSeq_iv(y) for y in transcript_intron_coords]

    flag=False
    for interval in transcript_intron_ivs:
        if interval.end-interval.start<1:
            print "transcript", transcript," has a problem with one or more of its introns:"
            print transcript_intron_coords
            flag = True
    if flag==False: introns_di[transcript.getName()]=transcript_intron_ivs

t2=time.time()
print "Execution time=", int(t2-t1), " seconds"

print len(transcript_feature_di), " transcript intervals defined"
print len(exons_di)," exon feature intervals defined"
print len(introns_di)," intron feature intervals defined"

df_total=pd.DataFrame()
df_cleaved_exon_counts=pd.DataFrame()
df_intron_counts=pd.DataFrame()
df_cleaved_intron_counts=pd.DataFrame()
df_CDS_total=pd.DataFrame()

column_headers=[]
p=multiprocessing.Pool(processes=24)

params_total_counts=[]
params_cleaved_exon_counts=[]
params_intron_counts=[]
params_cleaved_intron_counts=[]
params_cleaved_intron_counts=[]
params_CDS_counts=[]

for bamfilereader in BAM_file_readers:
    column_headers.append(bamfilereader.filename.split("/")[-1])
    params_total_counts.append((bamfilereader,transcript_feature_di))
    params_cleaved_exon_counts.append((bamfilereader,exons_di))
    params_intron_counts.append((bamfilereader,introns_di))
    params_cleaved_intron_counts.append((bamfilereader,introns_di))
    params_CDS_counts.append((bamfilereader,CDS_feature_di))

total_results=p.map(gobsChr.count_reads,params_total_counts)

print len(total_results)," Total Results Produced"

cleaved_exon_count_results=p.map(gobsChr.count_splicing_intermediates,params_cleaved_exon_counts)

print len(cleaved_exon_count_results), "Cleaved Exon Count Results Produced"

intron_count_results=p.map(gobsChr.count_reads_per_region,params_intron_counts)

print len(intron_count_results), "Intron Count Results Produced"

cleaved_intron_count_results=p.map(gobsChr.count_splicing_intermediates,params_cleaved_intron_counts)

print len(cleaved_intron_count_results), "Cleaved Intron Count Results Produced"

CDS_results=p.map(gobsChr.count_reads,params_CDS_counts)

print len(CDS_results), "CDS Count Results Produced"

for result, column_header in zip(total_results,column_headers):
    ds=pd.Series(result)
    df_total[column_header]=ds

for result,column_header in zip(cleaved_exon_count_results,column_headers):
    index=pd.MultiIndex.from_tuples(result.keys(),names=["Transcript","Intron"])
    ds=pd.Series(result,index=index)
    df_cleaved_exon_counts[column_header]=ds

for result, column_header in zip(intron_count_results, column_headers):
    index=pd.MultiIndex.from_tuples(result.keys(),names=["Transcript","Intron"])
    ds=pd.Series(result,index=index)
    df_intron_counts[column_header]=ds

for result, column_header in zip(cleaved_intron_count_results, column_headers):
    index=pd.MultiIndex.from_tuples(result.keys(),names=["Transcript","Intron"])
    ds=pd.Series(result,index=index)
    df_cleaved_intron_counts[column_header]=ds

for result, column_header in zip(CDS_results, column_headers):
    ds=pd.Series(result)
    df_CDS_total[column_header]=ds

df_total_sorted=df_total.sort()

df_cleaved_exon_counts_sorted=df_cleaved_exon_counts.sort()

df_intron_counts_sorted=df_intron_counts.sort()

df_cleaved_intron_counts_sorted=df_cleaved_intron_counts.sort()

df_CDS_total_sorted=df_cleaved_intron_counts.sort()


df_total_sorted.to_csv("{0}_totalcounts.tsv".format(sys.argv[2]),sep="\t")

df_cleaved_exon_counts_sorted.to_csv("{0}_cleaved_five_prime_exon_counts.tsv".format(sys.argv[2]),sep="\t")

df_intron_counts_sorted.to_csv("{0}_intron_counts.tsv".format(sys.argv[2]),sep="\t")

df_cleaved_intron_counts_sorted.to_csv("{0}_cleaved_intron_counts.tsv".format(sys.argv[2]),sep="\t")

df_CDS_total_sorted.to_csv("{0}_CDS_totals.tsv".format(sys.argv[2]),sep="\t")