'''Usage: python process_juncbase.py juncbase_output gff3_file fasta_file prefix
Please note: uses chr# format for chromosome names.'''

import sys
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome')
import SPTools as SP

juncbase_output = sys.argv[1]
gff3_file = sys.argv[2]
print gff3_file
fasta_file = sys.argv[3]
print fasta_file
prefix = sys.argv[4]

fasta_dict =  SP.make_fasta_dict(fasta_file)
junc_df, sample_list = SP.read_juncbase_output(juncbase_output)
seq_df = SP.get_junction_sequence(junc_df, gff3_file, fasta_dict)
pos_matrix_5prime, pos_matrix_3prime = SP.generate_consensus_matrix(gff3_file, fasta_dict, PSSM=True)
scored_df = SP.score_new_sites(seq_df, pos_matrix_5prime, pos_matrix_3prime, PSSM=True)
filt_df1, filt_df2 = SP.reformat_df(scored_df, sample_list)

#intron_ret_df = filt_df[filt_df['as_event_type'] == 'intron_retention']
#intron_ret_df = intron_ret_df.reset_index(drop=True)
#alt_donor = filt_df[filt_df['as_event_type'] == 'alternative_donor']
#alt_acceptor = filt_df[filt_df['as_event_type'] == 'alternative_acceptor']

filt_df1.to_csv('{0}_seq_score.tsv'.format(prefix), sep='\t', float_format='%.2f')
filt_df2.to_csv('{0}_PSI.tsv'.format(prefix), sep='\t', float_format='%.2f')
