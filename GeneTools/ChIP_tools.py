import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.append('/home/jordan/RNA-is-awesome/')
import GeneTools as GT
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/CHIP/')
sys.path.append('/home/jordan/RNA-is-awesome/CHIP/')
import Compare_RPKM
import pandas as pd
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt

def make_name(string):
    name = string.split('/')[-1].split('_sorted')[0]
    return name

def prep_bam(df, bam, tx_dict):
    name = bam.split('/')[-1].split('_sorted')[0]
    read_list = Compare_RPKM.count_reads_in_ChIP(tx_dict, bam)
    df[name] = read_list
    return df

def plot_min_max(lists, ax):
    all_min = 0.9*min([x for y in lists for x in y])
    all_max = 1.1*max([x for y in lists for x in y])
    ax.set_ylim([all_min,all_max])
    ax.set_xlim([all_min,all_max])
    ax.plot([all_min,all_max],[all_min,all_max],color='black')
    return ax

def ChIP_rpkm_scatter(WCE1_bam, WCE2_bam, WT1_bam, WT2_bam, Mut1_bam, Mut2_bam, gff3, Z_change=False, cen_tel=False):
    tx_dict = GT.build_transcript_dict(gff3)
    tx_dict = Compare_RPKM.make_promoter_dict(tx_dict, '/home/jordan/GENOMES/H99_chrom_lengths.json')
    
    df = pd.DataFrame(index=tx_dict.keys())
    bam_list = [WCE1_bam, WCE2_bam, WT1_bam, WT2_bam, Mut1_bam, Mut2_bam]
    for bam in bam_list:
        df = prep_bam(df, bam, tx_dict)
    
    for column in df:
        df[column] = pd.to_numeric(df[column])
    
    names = df.columns
    df[names[2]+' Normalized'] = df[names[2]]/df[names[0]]
    df[names[3]+' Normalized'] = df[names[3]]/df[names[1]]
    df[names[4]+' Normalized'] = df[names[4]]/df[names[0]]
    df[names[5]+' Normalized'] = df[names[5]]/df[names[1]]
    
    df['Enrichment 1'] = df[names[4]+' Normalized']/df[names[2]+' Normalized']
    df['Enrichment 2'] = df[names[5]+' Normalized']/df[names[3]+' Normalized']
    
    df = df[(df[names[2]] > 0) & (df[names[3]] > 0)]
    
    for_plot = []
    Z_scores = []
    for column in df:
        if 'Normalized' in column:
            print column
            df['log2 '+column.split(' ')[0]] = df[column].apply(np.log2)
            for_plot.append('log2 '+column.split(' ')[0])
        elif 'Enrichment' in column:
            df['Z-score '+column] = pd.Series(stats.mstats.zscore(df[column]), index=df.index)
            Z_scores.append('Z-score '+column)
    
    ## make the plot
    f, ax = plt.subplots(2, 2, figsize=(10,10))
    for axis in ax:
        axis[0] = plot_min_max([df[for_plot[0]].tolist()+df[for_plot[1]].tolist()+df[for_plot[2]].tolist()+df[for_plot[3]].tolist()], axis[0])
        axis[1] = plot_min_max([df[for_plot[0]].tolist()+df[for_plot[1]].tolist()+df[for_plot[2]].tolist()+df[for_plot[3]].tolist()], axis[1])
    
    ax[0,0].set_xlabel(for_plot[0])
    ax[0,0].set_ylabel(for_plot[1])
    ax[0,0].plot(df[for_plot[0]],df[for_plot[1]],'o', alpha=0.5, color='0.4')
    
    ax[0,1].set_xlabel(for_plot[2])
    ax[0,1].set_ylabel(for_plot[3])
    ax[0,1].plot(df[for_plot[2]],df[for_plot[3]],'o', alpha=0.5, color='0.4')
    
    ax[1,0].set_xlabel(for_plot[0])
    ax[1,0].set_ylabel(for_plot[2])
    
    if Z_change is True:
        ax[1,0].plot(df[for_plot[0]],df[for_plot[2]],'o', alpha=0.5, color='0.4')
        dep_df1 = df[(df[Z_scores[0]] >= 1.8) | (df[Z_scores[0]] <= -1.8)]
        ax[1,0].plot(dep_df1[for_plot[0]],dep_df1[for_plot[2]],'o', alpha=0.6, color='darkorchid')
    
    ax[1,1].set_xlabel(for_plot[1])
    ax[1,1].set_ylabel(for_plot[3])
    
    if Z_change is True:
        ax[1,1].plot(df[for_plot[1]],df[for_plot[3]],'o', alpha=0.5, color='0.4')
        dep_df2 = df[(df[Z_scores[1]] >= 1.8) | (df[Z_scores[1]] <= -1.8)]
        ax[1,1].plot(dep_df2[for_plot[1]],dep_df2[for_plot[3]],'o', alpha=0.6, color='darkorchid')
        
    if cen_tel is True:
        cen_df = df[df.index.str.contains('Cen')]
        tel_df = df[df.index.str.contains('tel')]
        ax[1,0].plot(cen_df[for_plot[0]], cen_df[for_plot[2]], 'o', alpha=0.8, color='crimson')
        ax[1,0].plot(tel_df[for_plot[0]], tel_df[for_plot[2]], 'o', alpha=0.8, color='mediumblue')
        
        ax[1,1].plot(cen_df[for_plot[1]], cen_df[for_plot[3]], 'o', alpha=0.8, color='crimson')
        ax[1,1].plot(tel_df[for_plot[1]], tel_df[for_plot[3]], 'o', alpha=0.8, color='mediumblue')

    plt.show()
    plt.clf()
    
