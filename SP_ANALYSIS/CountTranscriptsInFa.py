__author__ = 'jordanburke'

import sys
import collections
import operator

fin = open(sys.argv[1], "r")

fout = open("{0}_tally.fasta".format(sys.argv[1].split(".")[0]), "w")

transcript_list = []
for j, line in enumerate(fin):
    if j % 2 == 0:
        transcript_list.append(str(line.split(";")[0]).strip(">"))
    print transcript_list

counter = collections.Counter(transcript_list)
counter.most_common()
print counter

final_transcript_list = []
transcript_count = []
for transcript, count in counter.iteritems():
    final_transcript_list.append(transcript)
    transcript_count.append(count)

while n < len(transcript_count):
    line_list = [final_transcript_list[n],"\t",transcript_count[n],"\n"]
    line = "".join(line_list)
    fout.write(line)

fin.close()
fout.close()







