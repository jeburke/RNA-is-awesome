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

def ChIP_rpkm_scatter(WCE_bam, WT1_bam, WT2_bam, Mut1_bam, Mut2_bam, gff3, plot_name, Z_change=False, cen_tel=False):
    tx_dict = GT.build_transcript_dict(gff3)
    tx_dict = Compare_RPKM.make_promoter_dict(tx_dict, '/home/jordan/GENOMES/H99_chrom_lengths.json')
    
    df = pd.DataFrame(index=tx_dict.keys())
    bam_list = [WCE_bam, WT1_bam, WT2_bam, Mut1_bam, Mut2_bam]
    
    for bam in bam_list:
        df = prep_bam(df, bam, tx_dict)
    
    for column in df:
        df[column] = pd.to_numeric(df[column])
    
    names = df.columns
    df[names[1]+' Normalized'] = df[names[1]]/df[names[0]]
    df[names[2]+' Normalized'] = df[names[2]]/df[names[0]]
    df[names[3]+' Normalized'] = df[names[3]]/df[names[0]]
    df[names[4]+' Normalized'] = df[names[4]]/df[names[0]]
    
    df['Enrichment 1'] = df[names[3]+' Normalized']/df[names[1]+' Normalized']
    df['Enrichment 2'] = df[names[4]+' Normalized']/df[names[2]+' Normalized']
    
    df = df[(df[names[2]] > 0) & (df[names[3]] > 0)]
    
    for_plot = []
    for_plot2 = []
    Z_scores = []
    for column in df:
        if 'Normalized' in column:
            for_plot2.append(column)
            df['log2 RPKM '+column.split(' ')[0]] = df[column].apply(np.log2)
            for_plot.append('log2 RPKM '+column.split(' ')[0])
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
    ax[0,0].plot(df[for_plot[0]],df[for_plot[1]],'o', alpha=0.8, color='0.4')
    
    ax[0,1].set_xlabel(for_plot[2])
    ax[0,1].set_ylabel(for_plot[3])
    ax[0,1].plot(df[for_plot[2]],df[for_plot[3]],'o', alpha=0.8, color='0.4')
    
    ax[1,0].set_xlabel(for_plot[0])
    ax[1,0].set_ylabel(for_plot[2])
    
    if Z_change is True:
        ax[1,0].plot(df[for_plot[0]],df[for_plot[2]],'o', alpha=0.5, color='0.4')
        dep_df1 = df[(df[Z_scores[0]] >= 2) | (df[Z_scores[0]] <= -2)]
        ax[1,0].plot(dep_df1[for_plot[0]],dep_df1[for_plot[2]],'o', alpha=0.6, color='darkorchid')
    
    ax[1,1].set_xlabel(for_plot[1])
    ax[1,1].set_ylabel(for_plot[3])
    
    if Z_change is True:
        ax[1,1].plot(df[for_plot[1]],df[for_plot[3]],'o', alpha=0.5, color='0.4')
        dep_df2 = df[(df[Z_scores[1]] >= 2) | (df[Z_scores[1]] <= -2)]
        ax[1,1].plot(dep_df2[for_plot[1]],dep_df2[for_plot[3]],'o', alpha=0.6, color='darkorchid')
    
        up = set(df.index)
        down = set(df.index)
        for Z_col in Z_scores:
            up = up.intersection(df[df[Z_col] >= 2].index)
            down = down.intersection(df[df[Z_col] <= -2].index)
        df['Z-change'] = 'Unchanged'
        df.loc[up, 'Z-change'] = 'Up'
        df.loc[down, 'Z-change'] = 'Down'
    
    df.to_csv(plot_name+'.csv')    

    if cen_tel is True:
        cen_df = df[df.index.str.contains('Cen')]
        tel_df = df[df.index.str.contains('tel')]
        ax[1,0].plot(cen_df[for_plot[0]], cen_df[for_plot[2]], 'o', alpha=0.8, color='crimson', label='Centromeres')
        ax[1,0].plot(tel_df[for_plot[0]], tel_df[for_plot[2]], 'o', alpha=0.8, color='mediumblue', label='Telomeres')
        ax[1,0].legend()
        
        ax[1,1].plot(cen_df[for_plot[1]], cen_df[for_plot[3]], 'o', alpha=0.8, color='crimson', label='Centromeres')
        ax[1,1].plot(tel_df[for_plot[1]], tel_df[for_plot[3]], 'o', alpha=0.8, color='mediumblue', label='Telomeres')
        ax[1,1].legend()

    f.savefig(plot_name+'.eps', format='eps')
    print "Saved as "+plot_name+'.eps'
    
    plt.show()
    plt.clf()
    
    if cen_tel is True:
        cen_df = cen_df.sort_index()
        tel_df = tel_df.sort_index()
        
        cen_df['WT avg'] = (cen_df[for_plot2[0]] + cen_df[for_plot2[1]]).divide(2)
        cen_df['WT range'] = (cen_df[for_plot2[0]] - cen_df[for_plot2[1]]).apply(abs).divide(2)
        cen_df['Mut avg'] = (cen_df[for_plot2[2]] + cen_df[for_plot2[3]]).divide(2)
        cen_df['Mut range'] = (cen_df[for_plot2[2]] - cen_df[for_plot2[3]]).apply(abs).divide(2)
        
        tel_df['WT avg'] = (tel_df[for_plot2[0]] + tel_df[for_plot2[1]]).divide(2)
        tel_df['WT range'] = (tel_df[for_plot2[0]] - tel_df[for_plot2[1]]).apply(abs).divide(2)
        tel_df['Mut avg'] = (tel_df[for_plot2[2]] + tel_df[for_plot2[3]]).divide(2)
        tel_df['Mut range'] = (tel_df[for_plot2[2]] - tel_df[for_plot2[3]]).apply(abs).divide(2)
        
        f2, ax2 = plt.subplots(2, 1, figsize=(10,10))
        width = 0.35
        
        ind_cen = np.arange(len(cen_df))
        ax2[0].bar(ind_cen, cen_df['WT avg'], width, color='mediumblue', yerr=cen_df['WT range'], label="WT")
        ax2[0].bar(ind_cen + width, cen_df['Mut avg'], width, color='crimson', yerr=cen_df['Mut range'], label="Mutant")
        
        ax2[0].set_ylabel('RPKM')
        ax2[0].set_title('Centromeres')
        ax2[0].set_xticks(ind_cen + width / 2)
        ax2[0].set_xticklabels(cen_df.index)
        ax2[0].legend()
        
        ind_tel = np.arange(len(tel_df))
        ax2[1].bar(ind_tel, tel_df['WT avg'], width, color='mediumblue', yerr=tel_df['WT range'], label="WT" )
        ax2[1].bar(ind_tel + width, tel_df['Mut avg'], width, color='crimson', yerr=tel_df['Mut range'], label="Mutant")
        
        ax2[1].set_ylabel('RPKM')
        ax2[1].set_title('Telomeres')
        ax2[1].set_xticks(ind_tel + width / 2)
        ax2[1].set_xticklabels(tel_df.index, rotation='vertical')
        ax2[1].legend()

        f2.savefig(plot_name+'_bar.eps', format='eps')
        print "Saved as "+plot_name+'_bar.eps'

        plt.show()
        plt.clf()
