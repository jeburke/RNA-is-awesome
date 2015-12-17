__author__ = 'jordanburke'

import sys
import collections
import operator

fin = open(sys.argv[1], "r")

fout = open("{0}_tally.fasta".format(sys.argv[1].split(".")[0]), "w")

transcript_list = []
for line in fin:
    if line.startswith(">"):
        transcript_list.append(str(line.split(";")[0]).strip(">"))
print transcript_list[20]

counter = collections.Counter(transcript_list)
counter_sorted = counter.most_common()
#print counter

final_transcript_list = []
transcript_count = []

a=0
while a < len(counter_sorted):
    final_transcript_list.append(counter_sorted[a][0])
    transcript_count.append(counter_sorted[a][1])
    a += 1

b = 0
while b < len(transcript_count):
    line_list = [str(final_transcript_list[b]),"\t",str(transcript_count[b]),"\n"]
    line = "".join(line_list)
    fout.write(line)
    b += 1

fin.close()
fout.close()







