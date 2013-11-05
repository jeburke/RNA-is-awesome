import numpy
import sys
#import pylab

genome = {}

#Digest the genome into fragments using 1 or 2 restriction sites (of any length)

def digest(genome, site1, site2=None):
        fragments = []
        fragsize = []
        for chrom, sequence in genome.iteritems():
                fragments = sequence.split(site1)
                if site2 == None:
                        for x in fragments:
                                fragsize.append(len(x))
                else:
                        for x in fragments:
                                for y in x.split(site2):
                                        fragsize.append(len(y))
                return fragsize

#Compile histogram of fragment sizes

def histogram(fragsize):
        bins = 20
        hist, bin_edges = numpy.histogram(fragsize, bins=bins)
        # hist: [number in bin]
        # bin_edges: [edges of bin]
        return hist, bin_edges

#       P.figure()
#       hist, bin_edges, patches = P.hist(x, bin_edges, normed=1, histtype = 'bar', rwidth=0.8)
#       P.show()

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
for key, value in genome.iteritems():
        print key
        print len(str(value))
print str(len(genome))

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

fout.close()
