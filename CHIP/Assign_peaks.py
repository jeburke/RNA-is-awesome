import sys
import json
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneTools as GT
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/SP_ANALYSIS')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/SP_ANALYSIS')
import SPTools as SP
import subprocess
import shlex


def sort_tx_by_chrom(tx_dict):
    transcript_by_chr = {}
    for transcript, coords in tx_dict.iteritems():
        if transcript[-2] == 'T':
            transcript = transcript[:-2]
        chromosome = coords[3]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(transcript)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(transcript)
    return transcript_by_chr

def make_intergenic_dict(tx_dict, tx_by_chrom):
    print "Generating intergenic dictionary"
    int_dict = {}
    for chrom, tx_list in tx_by_chrom.iteritems():
        tx_list = [x for x in tx_list if x.endswith('T0')]
        chrom_df = pd.DataFrame(index=range(len(tx_list)), columns=['transcript','start','end','strand'])
        starts = []
        ends = []
        strands = []
        for tx in tx_list:
            starts.append(tx_dict[tx][0])
            ends.append(tx_dict[tx][1])
            strands.append(tx_dict[tx][2])
        chrom_df['transcript'] = tx_list
        chrom_df['start'] = starts
        chrom_df['end'] = ends
        chrom_df['strand'] = strands
        
        chrom_df = chrom_df.sort_values('start')
        chrom_df = chrom_df.reset_index(drop=True)
        n=0
        for n in range(len(chrom_df)):
            if chrom_df.ix[n]['strand'] == '+' and n > 0:
                prev_end = chrom_df.ix[n-1]['end']
                if prev_end < chrom_df.ix[n]['start']:
                    int_dict[chrom_df.ix[n]['transcript']] = [chrom_df.ix[n-1]['end'],chrom_df.ix[n]['end'],chrom_df.ix[n]['strand'], chrom]
            elif chrom_df.ix[n]['strand'] == '-' and n < len(chrom_df)-1:
                next_start = chrom_df.ix[n+1]['start']
                if next_start > chrom_df.ix[n]['end']:
                    int_dict[chrom_df.ix[n]['transcript']] = [chrom_df.ix[n]['start'],chrom_df.ix[n+1]['start'],chrom_df.ix[n]['strand'], chrom]
                    
    new_tx_dict = {k:v for k, v in tx_dict.items() if k not in int_dict}
    int_dict.update(new_tx_dict)
    return int_dict

def make_promoter_dict(tx_dict, chrom_lengths):
    if type(chrom_lengths) == str:
        with open(chrom_lengths, 'r') as f:
            chrom_lengths = json.load(f)
    prom_dict = {}
    for tx, info in tx_dict.iteritems():
        if tx[-2] == 'T':
            tx = tx[:-2] 
        if info[2] == '+':
            if info[0] >= 1000:
                new_start = info[0]-1000
            else:
                new_start = 0
            prom_dict[tx] = [new_start, info[1], info[2], info[3]]
        elif info[2] == '-':
            if info[1] <= chrom_lengths[info[3]]-1000:
                new_end = info[1]+1000
            else:
                new_end = chrom_lengths[info[3]]
            prom_dict[tx] = [info[0], new_end, info[2], info[3]]
    return prom_dict
            

def assign_peak_to_tx(tx_by_chrom, tx_dict, MACS_out, cutoff=None):
    print "Assigning peaks to genes"
    df = pd.read_csv(MACS_out, sep='\t', skiprows=24)
    transcripts = []
    for index, row in df.iterrows():
        row_tx = ''
        for tx in tx_by_chrom[row['chr']]:
            if row['abs_summit'] >= tx_dict[tx][0] and row['abs_summit'] <= tx_dict[tx][1]:
                if row_tx != '':
                    row_tx = row_tx+','+tx
                else:
                    row_tx = tx
        transcripts.append(row_tx)
    df['transcript'] = transcripts
    df.index = df['transcript']
    if cutoff is not None:
        df = df[df['fold_enrichment'] > cutoff]
    print len(df)
    return df

def find_best_peaks(df, tx_dict, max_genes=None):
    new_df = pd.DataFrame(index = tx_dict.keys(), columns = df.columns)
    for tx in tx_dict:
        tx_df = df[df.index.str.contains(tx)]
        if len(tx_df > 0):
            tx_df = tx_df.sort_values('abs_summit')
            tx_df = tx_df.reset_index(drop=True)
            best_row = tx_df.iloc[0].tolist()
            new_df.loc[tx] = best_row
    new_df = new_df.dropna()
    if max_genes is not None:
        new_df = new_df.sort_values('-log10(pvalue)', ascending=False)
        new_df = df.iloc[:max_genes]
    return new_df

def generate_sequence_file(df, tx_dict, fasta, prefix):
    with open(fasta,'r') as f:
        fasta = json.load(f)
    seq_dict = {}
    for tx, r in df.iterrows():
        for each_tx in tx.split(','):
            if len(each_tx) > 1:
                name = each_tx+':'+str(r['abs_summit'])
                seq = GT.seq_simple(r['chr'], r['start'], r['end'], tx_dict[each_tx][2], fasta)
                seq_dict[name] = seq
        
    with open(prefix+'_peak_sequences.fa', 'w') as fout:
        for name, seq in seq_dict.iteritems():
            fout.write('>'+name+'\n')
            fout.write(seq+'\n')
            
def split_by_gene(df, gene_list_file):
    gene_list = []
    with open(gene_list_file,'r') as f:
        for line in f:
            gene_list.append(line.strip())
    
    in_list_df = df[df.index.isin(gene_list)]
    other_df = df[~df.index.isin(gene_list)]
    
    return in_list_df, other_df

def call_meme(prefix, minsites, split=False):
    if split is False:
        command = '/usr/local/bin/meme-chip -meme-maxw 8 -meme-minsites '+str(minsites)+' -oc '+prefix+' '+prefix+'_peak_sequences.fa'
        args = shlex.split(command)
        print command
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    else:
        command1 = '/usr/local/bin/meme-chip -noecho -meme-maxw 8 -meme-minsites '+str(minsites[0])+' -oc '+prefix+'_in_list '+prefix+'_in_list_peak_sequences.fa'
        args1 = shlex.split(command1)
        print command1
        p1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        command2 = '/usr/local/bin/meme-chip -noecho -meme-maxw 8 -meme-minsites '+str(minsites[1])+' -oc '+prefix+'_other '+prefix+'_other_peak_sequences.fa'
        args2 = shlex.split(command2)
        print command2
        p2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    

def main():
    gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    fasta = '/home/jordan/GENOMES/H99_fa.json'
    chrom_lengths = '/home/jordan/GENOMES/H99_chrom_lengths.json'
    prefix = sys.argv[1].split('/')[-1].split('.')[0]
    print prefix
    tx_dict = SP.build_transcript_dict(gff3)
    tx_by_chrom = sort_tx_by_chrom(tx_dict)
    int_dict = make_promoter_dict(tx_dict, chrom_lengths)
    peak_df = assign_peak_to_tx(tx_by_chrom, int_dict, sys.argv[1], cutoff=2)
    #peak_df = assign_peak_to_tx(tx_by_chrom, int_dict, sys.argv[1])
    peak_df = find_best_peaks(peak_df, int_dict, max_genes=300)
    if len(sys.argv) == 3:
        gene_list_file = sys.argv[2]
        in_list, other = split_by_gene(peak_df, gene_list_file)
        in_list.to_csv(prefix+'_by_gene_in_list.csv')
        other.to_csv(prefix+'_by_gene_other.csv')
        generate_sequence_file(in_list, int_dict, fasta, prefix+'_in_list')
        generate_sequence_file(other, int_dict, fasta, prefix+'_other')
        split = True
        minsites = [int(0.75*len(in_list)),int(0.75*len(other))]
        if minsites[0] > 600: minsites[0] = 600
        if minsites[1] > 600: minsites[1] = 600
    else:
        peak_df.to_csv(prefix+'_by_gene.csv')
        generate_sequence_file(peak_df, int_dict, fasta, prefix)
        split = False
        minsites = int(0.75*len(peak_df))
        if minsites > 600: minsites = 600
    call_meme(prefix, minsites, split=split)
    
if __name__ == "__main__":
    main()
        