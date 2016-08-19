'''Usage: python procues_juncbase.py juncbase_output gff3_file fasta_file prefix
Please note: uses chr# format for chromosomes.'''

import sys
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/SP_ANALYSIS/')
import SPTools as SP

juncbase_output = sys.argv[1]
gff3_file = sys.argv[2]
fasta_file = sys.argv[3]
prefix = sys.arv[4]

junc_df = SP.read_juncbase_output(juncbase_output)
seq_df = SP.get_junction_sequence(junc_df, gff3_file, fasta_file)
pos_matrix_5prime, pos_matrix_3prime = SP.generate_consensus_matrix()
scored_df = SP.score_new_sites(seq_df, pos_matrix_5prime, pos_matrix_3prime)
filt_df = SP.reformat_df(scored_df)

#intron_ret_df = filt_df[filt_df['as_event_type'] == 'intron_retention']
#intron_ret_df = intron_ret_df.reset_index(drop=True)
#alt_donor = filt_df[filt_df['as_event_type'] == 'alternative_donor']
#alt_acceptor = filt_df[filt_df['as_event_type'] == 'alternative_acceptor']

filt_df.to_csv('{0}.tsv'.format(prefix), sep='\t')