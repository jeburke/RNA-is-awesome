__author__ = 'jordanburke'

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import gobsChr
import pandas as pd
import time
import sys

t1=time.time()
basedir="./"
print 'reading configuration file'

#############################################
#  This code reads a configuration file     #
#  that contains the paths to the .bam files#
#  to be analyzed and the GFF3 file used to #
#  generate those bamfiles using TOPHAT     #
#############################################

config_file=open(basedir+"conf/LaSSO.config")
for line in config_file:
    data=line.split("=")
    if data[0]=="GFF":
        GFF_file_name=data[1].strip()
        print GFF_file_name

transcripts=gobsChr.transcript_factory(basedir+"conf/H99.config")


transcript_feature_di={}
transcript_len_di={}
introns_di={}

for transcript in transcripts:
    transcript_coord=transcript.getCoordinates()
    transcript_iv=gobsChr.HTSeq_iv(transcript_coord)
    transcript_feature_di[transcript.getName()]=transcript_iv
    transcript_len_di[transcript.getName()]=transcript_iv.end-transcript_iv.start
    transcript_intron_coords = transcript.getAllIntronCoordinates()
    if len(transcript_intron_coords)==0:
        print "no introns in ",transcript
        continue
    if transcript_intron_coords[0].getStrand() == "Watson":
        transcript_intron_coords.sort()
    else:
        transcript_intron_coords.sort(reverse=True)
    transcript_intron_ivs=[gobsChr.HTSeq_iv(x) for x in transcript_intron_coords]
    flag=False
    for interval in transcript_intron_ivs:
        if interval.end-interval.start<1:
            print "transcript", transcript," has a probem with one or more of its exons:"
            print transcript_intron_coords
            flag=True
    if flag==False: introns_di[transcript.getName()]=transcript_intron_ivs
#    print transcript.getName()
t2=time.time()

print len(transcript_feature_di), " transcript intervals defined"
print len(introns_di)," intron feature intervals defined"


transcript_list = []
chrom_list = []
start_list = []
end_list = []
strand_list = []
exon_number = []

print introns_di["CNAG_04612T2"][0].chrom
print introns_di["CNAG_04612T0"][1].start
print introns_di["CNAG_04612T0"][1].end
print introns_di["CNAG_04612T0"][1].strand

#print introns_di["CNAG_07977TO"][0].start

for gene, intron_list in introns_di.iteritems():
    n = 0
    while n < len(intron_list):
        transcript_list.append(gene)
        chrom_list.append(introns_di[gene][n].chrom)
        start_list.append(introns_di[gene][n].start)
        end_list.append(introns_di[gene][n].end)
        strand_list.append(introns_di[gene][n].strand)
        n+=1
        exon_number.append(n)

strand_number = []
for strand in strand_list:
    if strand == "+":
        strand_number.append(1)
    elif strand == "-":
        strand_number.append(-1)

print len(transcript_list)
print len(chrom_list)
print len(start_list)
print len(end_list)
print len(strand_list)

print transcript_list[1]
print exon_number[1]
print exon_number[10]
print transcript_list[10]
print transcript_list[20]

print start_list[1]
print start_list[20]
print end_list[20]
print strand_list[20]
print strand_number[20]
print exon_number[20]

fout = open("{0}_LaSSO_header.txt".format(sys.argv[1]), "w")
a = 0
while a < len(transcript_list):
    line_list = [str(transcript_list[a]),";",str(chrom_list[a]),";",str(strand_number[a]),";",str(start_list[a]),";",str(end_list[a]),";",str(transcript_list[a]),":exon:",str(exon_number[a]),"-",str(transcript_list[a]),":exon:",str(exon_number[a]+1),"\n"]
    line = "".join(line_list)
    fout.write(line)
    a += 1
fout.close()