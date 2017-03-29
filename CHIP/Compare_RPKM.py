import sys
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneTools as GT
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import json
import pysam

def make_promoter_dict(tx_dict, chrom_lengths):
    if type(chrom_lengths) == str:
        with open(chrom_lengths, 'r') as f:
            chrom_lengths = json.load(f)
    prom_dict = {}
    for tx, info in tx_dict.iteritems():
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
    prom_dict = {k[:-2]:v for k, v in prom_dict.items() if k.endswith('T0')}
    return prom_dict

def count_reads_in_ChIP(prom_dict, bam_file):
    bam = pysam.Samfile(bam_file)
    total = GT.count_aligned_reads(bam_file)
    print total
    
    rpkm_s = pd.Series(index=prom_dict.keys())
    for tx, info in prom_dict.iteritems():
        read_count = 0
        if info[1]-info[0] <= 0:
            pass
        else:
            read_count += GT.count_reads_in_window(bam, info[3], info[0], info[1], '+')
            read_count += GT.count_reads_in_window(bam, info[3], info[0], info[1], '-')
            length = (info[1]-info[0])/1000.
            rpkm = read_count/length/total
            rpkm_s[tx] = rpkm
    rpkm_s = rpkm_s.dropna()
    return rpkm_s

def compare_reads(ChIP_bam1, WCE_bam1, ChIP_bam2, WCE_bam2, prom_dict, prefix, fit_reg=False):
    name1 = ChIP_bam1.split('/')[-1].split('_sorted')[0]
    print name1
    sample1 = count_reads_in_ChIP(prom_dict, ChIP_bam1)
    print sum(sample1)
             
    name2 = ChIP_bam2.split('/')[-1].split('_sorted')[0]
    print '\n'+name2
    sample2 = count_reads_in_ChIP(prom_dict, ChIP_bam2)
    print sum(sample2)
    
    WCE_name1 = WCE_bam1.split('/')[-1].split('_sorted')[0]
    print '\n'+WCE_name1
    WCE1 = count_reads_in_ChIP(prom_dict, WCE_bam1)
    print sum(WCE1)
    
    WCE_name2 = WCE_bam2.split('/')[-1].split('_sorted')[0]
    print '\n'+WCE_name2
    WCE2 = count_reads_in_ChIP(prom_dict, WCE_bam2)
    print sum(WCE2)
    
    df = pd.DataFrame(index=prom_dict.keys(), columns=[name1, name2, WCE_name1, WCE_name2, name1+' Normalized', name2+' Normalized', 'Enrichment'])
    
    df[name1] = sample1
    df[name2] = sample2
    df[WCE_name1] = WCE1
    df[WCE_name2] = WCE2
    df[name1+' Normalized'] = df[name1]/df[WCE_name1]
    df[name2+' Normalized'] = df[name2]/df[WCE_name2]
    df['Enrichment'] = df[name2+' Normalized']/df[name1+' Normalized']
    
    df = df[df[name2+' Normalized'] > 2]
    print len(df)
    
    fig = plt.figure(figsize=(8,6))
    ax = sns.regplot(x=name1+' Normalized', y=name2+' Normalized', data=df, fit_reg=fit_reg)
    plt.show()
    fig.save_fig(prefix+'_scatter.pdf',format='pdf')
    plt.clf()
    
    return df
        
def main():
    ChIP_bam1 = sys.argv[1]
    WCE_bam1 = sys.argv[2]
    ChIP_bam2 = sys.argv[3]
    WCE_bam2 = sys.argv[4]
    try:
        prefix = sys.argv[5]
    except IndexError:
        prefix = 'Sample1_vs_Sample2'
    gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    tx_dict = GT.build_transcript_dict(gff3)
    prom_dict = make_promoter_dict(tx_dict, '/home/jordan/GENOMES/H99_chrom_lengths.json')
    rpkm_df = compare_reads(ChIP_bam1, WCE_bam1, ChIP_bam2, WCE_bam2, prom_dict, prefix, fit_reg=True)
    rpkm_df.to_csv(prefix+'rpkm_comparison.csv')
    
if __name__ == "__main__":
    main()
        