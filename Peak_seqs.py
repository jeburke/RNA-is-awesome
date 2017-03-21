'''Usage: python Peak_seqs.py input_csv window_size
Input file should be excel csv saved as Windows csv
Columns are transcript, chromosome, peak center
window_size is the total window for the sequence (so window/2 on either side of peak center'''

import GeneTools as GT
import sys
import json

input_file = sys.argv[1]
with open('/home/jordan/GENOMES/H99_fa.json','r') as fa:
    fasta_file = json.load(fa)
gff3_file = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
window = int(sys.argv[2])
seq_list = GT.get_peak_sequence(input_file, fasta_file, gff3_file, window=window)
