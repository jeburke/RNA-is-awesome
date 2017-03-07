import GeneTools as GT
import sys

input_file = sys.argv[1]
fasta_file = '/home/jordan/GENOMES/CNA3-gobs.fa'
gff3_file = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
window = int(sys.argv[2])
seq_list = GT.get_peak_sequence(input_file, fasta_file, gff3_file, window=1000)