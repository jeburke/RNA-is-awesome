'''Created on 4 November 2013

Usage: python DigestGenome.py <genome>.fasta <prefix> <CUTSITE1> <CUTSITE2(optional)>
Sequence of restriciton site should be in all caps. 2nd restriction site is optional.
@author: jordan
'''

import numpy
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path

genome = {}

#Digest the genome into fragments using 1 or 2 restriction sites (of any length)

def digest(genome, site1, site2=None):
        fragments = []
        fragsize = []
	less50 = 0
        for chrom, sequence in genome.iteritems():
                fragments = sequence.split(site1)
                if site2 == None:
                        for x in fragments:
                                fragsize.append(len(x))
                else:
                        for x in fragments:
                                for y in x.split(site2):
                                        fragsize.append(len(y))
		for i in range(len(fragsize)):
			if fragsize[i] < 50:
				less50 += 1
		print "Number of fragments less than 50 bp: %d" % (less50)
                return fragsize

#Compile histogram of fragment sizes

def histogram(fragsize):
#	bins = 40
	bins = [0, 25, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400]
        hist, bin_edges = numpy.histogram(fragsize, bins=bins)
        # hist: [number in bin]
	# bin_edges: [edges of bin]
	
	n, b, patches = plt.hist(fragsize, bins, normed=1, facecolor='blue', alpha=0.5)
	plt.gca().set_xscale("log")
	plt.xlabel('Fragment sizes')
	plt.ylabel('Number')
	plt.subplots_adjust(left=0.15)	

	plt.savefig('{0}_genome_digest.png'.format(sys.argv[2]), format = 'png')

        return hist, bin_edges

#Open the genome fasta file and create a dictionary containing each chromosome

def get_chrom(file):
        chromNum = 0
        genome[chromNum] = []
        fin = open(file, "r")
        for line in fin:
                if line.startswith(">"):
                        genome[chromNum] = ''.join(genome[chromNum])
                        chromNum += 1
                        genome[chromNum] = []
                else:
                        genome[chromNum].append(line.strip())
        genome[chromNum] = ''.join(genome[chromNum])
        fin.close()
        del genome[0]
        return genome

genome = get_chrom(sys.argv[1])
#for key, value in genome.iteritems():
#        print key
#        print len(str(value))
print "Number of chromosomes: %d" % (len(genome))

if len(sys.argv) > 4:
        fragsize = digest(genome, sys.argv[3], sys.argv[4])
else:
        fragsize = digest(genome, sys.argv[3])

hist, bin_edges = histogram(fragsize)

fout = open("{0}_genome_digest.txt".format(sys.argv[2]), "w")
header_list = ["Fragment Size", "Number Fragments", "\n"]
header = "\t".join(header_list)
fout.write(header)
for i in range(len(hist)):
        line_list = [str(bin_edges[i+1]), str(hist[i]), "\n"]
        line = "\t".join(line_list)
        fout.write(line)

print "Histogram saved as: %s" %("{0}_genome_digest.txt".format(sys.argv[2]))
fout.close()
