'''Script to split transcripts in a fasta file either at the 5' or 3' end by some number of nucleotides
Usage: python Split_transcript_fasta.py <input fasta file> <five_prime or three_prime> <number of nucleotides>
Example: python Split_transcript_fasta.py W303_JRIU00000000_SGD_cds.fsa three_prime 100'''

import sys
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneTools
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/')
import SPTools as SP

cds_dict = GeneTools.read_fasta(sys.argv[1])
split_dict = GeneTools.split_cds(cds_dict, sys.argv[2], int(sys.argv[3]))
GeneTools.write_new_fasta(split_dict, sys.argv[1])