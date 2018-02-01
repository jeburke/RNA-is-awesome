import pysam
import pandas as pd
import numpy as np
import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
import SPTools as SP
from collections import OrderedDict
from matplotlib import pyplot as plt

def make_transcript_df(gff3):
    if 'pombe' in gff3.lower():
        organism='pombe'
    else:
        organism=None
    
    # Get transcript dictionary
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    
    # Organize by transcript
    tx_dict = OrderedDict(sorted(tx_dict.items(), key=lambda t: t[0]))
    
    # Convert to dataframe
    tx_df = pd.DataFrame(index=tx_dict.keys(), columns=['start','end','strand','chromosome'])
    for n, col in enumerate(tx_df.columns):
        tx_df.loc[:,col] = zip(*tx_dict.values())[n]
    
    # Add CDS starts and ends
    CDS_starts = [min(l) if len(l) > 0 else np.NaN for l in zip(*tx_dict.values())[4]]
    CDS_ends = [max(l) if len(l) > 0 else np.NaN for l in zip(*tx_dict.values())[5]]
    tx_df.loc[:,'CDS start'] = CDS_starts
    tx_df.loc[:,'CDS end'] = CDS_ends
    
    return tx_df

def count_siRNA_reads(bam_list, gff3):
    tx_df = make_transcript_df(gff3)

    for bam_file in bam_list:
        name = bam_file.split('/')[-1].split('_sorted.bam')[0]
        print name
        bam = pysam.Samfile(bam_file)
        sense = []
        antisense = []
        
        for ix, r in tx_df.iterrows():
            rc1 = 0
            rc2 = 0
            tx_reads = bam.fetch(r['chromosome'],r['start'],r['end'])
            for read in tx_reads:
                if r['strand'] == '+':
                    if not read.is_reverse:
                        rc1 += 1
                    elif read.is_reverse:
                        rc2 += 1
                elif r['strand'] == '-':
                    if read.is_reverse:
                        rc1 += 1
                    elif not read.is_reverse:
                        rc2 += 1
            sense.append(rc1)
            antisense.append(rc2)
        tx_df[name+' siRNA sense'] = sense
        tx_df[name+' siRNA antisense'] = antisense
        
    for column in tx_df.columns:
        if 'siRNA' in column:
            tx_df[column+' norm'] = tx_df[column].divide(sum(tx_df[column])/1000000.)
    
    # Filter based on a minimum number of reads in the antisense columns
    anti_cols = [x for x in tx_df.columns if x.endswith('antisense')]
    filtered = tx_df[anti_cols].sum(axis=1)
    filtered = filtered[filtered > 100*len(anti_cols)]
    tx_df = tx_df[tx_df.index.isin(filtered.index)]
    
    return tx_df

def rep_scatters(WTreps, MUTreps, df, name='Compare'):
    s1 = df[MUTreps[0]]/df[WTreps[0]]
    s2 = df[MUTreps[1]]/df[WTreps[1]]
    up = set(s1[s1>2].index).intersection(s2[s2>2].index)
    print "Up in mutant: "+str(len(up))
    with open(name+'_up.txt', 'w') as fout:
        for gene in up:
            fout.write(gene+'\n')
    with open(name+'_up_genes.txt', 'w') as fout:
        for gene in set([x[:-2] for x in up]):
            fout.write(gene+'\n')
    
    down = set(s1[s1<0.5].index).intersection(s2[s2<0.5].index)
    print "Down in mutant: "+str(len(down))
    with open(name+'_down.txt', 'w') as fout:
        for gene in down:
            fout.write(gene+'\n')
    with open(name+'_down_genes.txt', 'w') as fout:
        for gene in set([x[:-2] for x in down]):
            fout.write(gene+'\n')
    
    fig, ax = plt.subplots(ncols=2, figsize=(8,4), sharex=True, sharey=True)
    ax[0].scatter(df[WTreps[0]].apply(np.log2), df[MUTreps[0]].apply(np.log2), color='0.3', s=15)
    ax[1].scatter(df[WTreps[1]].apply(np.log2), df[MUTreps[1]].apply(np.log2), color='0.3', s=15)

    ax[0].scatter(df[(df.index.isin(up)) | (df.index.isin(down))][WTreps[0]].apply(np.log2), 
                  df[(df.index.isin(up)) | (df.index.isin(down))][MUTreps[0]].apply(np.log2), color='orchid', s=15)
    ax[1].scatter(df[(df.index.isin(up)) | (df.index.isin(down))][WTreps[1]].apply(np.log2), 
                  df[(df.index.isin(up)) | (df.index.isin(down))][MUTreps[1]].apply(np.log2), color='orchid', s=15)
    ax[0].set_xlabel(WTreps[0].split('siRNA')[0]+'RPM', fontsize=14)
    ax[0].set_ylabel(MUTreps[0].split('siRNA')[0]+'RPM', fontsize=14)
    ax[1].set_xlabel(WTreps[1].split('siRNA')[0]+'RPM', fontsize=14)
    ax[1].set_ylabel(MUTreps[1].split('siRNA')[0]+'RPM', fontsize=14)

    for sub in ax:
        sub.set_xlim(sub.get_xlim())
        sub.set_ylim(sub.get_xlim())
        sub.plot(sub.get_xlim(), sub.get_xlim(), '--', color='0.7', zorder=0)

    fig.tight_layout()
    plt.show()
    plt.clf()
    return fig, up, down

def read_size_distribution(df, bam_list, tx_list=None):
    if tx_list is not None:
        df = df[df.index.isin(tx_list)]
    
    read_sizes = {}
    for bam in bam_list:
        print bam
        read_sizes[bam] = []
        open_bam = pysam.Samfile(bam)
        for ix, r in df.iterrows():
            reads = open_bam.fetch(r['chromosome'],r['start'], r['end'])
            for read in reads:
                if read.is_reverse and r['strand'] == '+':
                    read_sizes[bam].append(len(read.query_alignment_sequence))
                elif not read.is_reverse and r['strand'] == '-':
                    read_sizes[bam].append(len(read.query_alignment_sequence))
    
    bins=np.arange(16,30)
    fig, ax = plt.subplots(len(bam_list), figsize=(4,2*len(bam_list)), sharex=True, sharey=True)
    for n, bam in enumerate(bam_list):
        ax[n].hist(read_sizes[bam], bins=bins, normed=True, label=bam.split('_sorted')[0].split('/')[-1], color='lavender', width=0.5)
        ax[n].legend()
    ax[0].set_xlim(16,30)
    
    fig.tight_layout()
    plt.show()
    plt.clf()
    return fig

def main():
    bam_list = []
    for arg in sys.argv[1:-1]:
        bam_list.append(arg)
    if len(bam_list)%4 != 0:
        "Uneven number of bam files provided!"
        return None
    
    bam_names = []
    for bam in bam_list:
        bam_names.append(bam.split('/')[-1].split('_sorted.bam')[0]+' siRNA antisense norm')
    
    name = sys.argv[-1]
    
    gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    df = count_siRNA_reads(bam_list, gff3)
    df.to_csv(name+'.csv')
    
    fig, up, down = rep_scatters(bam_names[0:2], bam_names[2:4], df, name=name)
    fig.savefig(name+'_A_scatters.pdf', format='pdf', bbox_inches='tight')
    fig = read_size_distribution(df, bam_list, tx_list=up)
    fig.savefig(name+'_A_up_read_dist.pdf', format='pdf', bbox_inches='tight')
    fig = read_size_distribution(df, bam_list, tx_list=down)
    fig.savefig(name+'_A_down_read_dist.pdf', format='pdf', bbox_inches='tight')
    
    if len(bam_list) >= 8:
        fig, up, down = rep_scatters(bam_names[4:6], bam_names[6:8], df, name=name)
        fig.savefig(name+'_B_scatters.pdf', format='pdf', bbox_inches='tight')
        fig = read_size_distribution(df, bam_list, tx_list=up)
        fig.savefig(name+'_B_up_read_dist.pdf', format='pdf', bbox_inches='tight')
        fig = read_size_distribution(df, bam_list, tx_list=down)
        fig.savefig(name+'_B_down_read_dist.pdf', format='pdf', bbox_inches='tight')
        
if __name__ == "__main__":
    main()