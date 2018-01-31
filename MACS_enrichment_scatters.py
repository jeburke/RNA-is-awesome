import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys

def prep_csv(csv):
    df = pd.read_csv(csv)
    df = df.sort_values(['fold_enrichment'], ascending=False)
    df = df.drop_duplicates(subset=['chr','start','end'])
    return df
    
def compare_peaks_MACS(df1, df2, min_overlap=0.5):
    enrich_mut1 = []
    enrich_mut2 = []
    for ix, r in df1.iterrows():
        match = False
        chrom_df = df2[df2['chr'] == r['chr']]
        range1 = set(range(r['start'],r['end']))
        for ix2, r2 in chrom_df.iterrows():
            range2 = range(r2['start'],r2['end'])
            if len(range1.intersection(range2))/float(len(range1)) >= min_overlap:
                enrich_mut1.append(r2['fold_enrichment'])
                enrich_mut2.append(r2['fold_enrichment2'])
                match = True
                break
        if match is False:
            enrich_mut1.append(np.NaN)
            enrich_mut2.append(np.NaN)
    df1.loc[:,'fold_enrichment mut1'] = enrich_mut1
    df1.loc[:,'fold_enrichment mut2'] = enrich_mut2
    return df1

def scatter_plots(df, csv1, csv2):
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10,10))
    ax[0][0].scatter(df['fold_enrichment'], df['fold_enrichment2'], s=15, color='0.3', alpha=0.5)
    ax[0][1].scatter(df['fold_enrichment mut1'], df['fold_enrichment mut2'], s=15, color='0.3', alpha=0.5)
    ax[1][0].scatter(df['fold_enrichment'], df['fold_enrichment mut1'], s=15, color='0.3', alpha=0.5)
    ax[1][1].scatter(df['fold_enrichment2'], df['fold_enrichment mut2'], s=15, color='0.3', alpha=0.5)

    min_xy = min(ax[0][0].get_xlim()[0], ax[0][0].get_ylim()[0])
    max_xy = max(ax[0][0].get_xlim()[1], ax[0][0].get_ylim()[1])
    
    for subax in ax:
        for sub in subax:
            min_xy = min(min_xy, sub.get_xlim()[0], sub.get_ylim()[0])
            max_xy = max(max_xy, sub.get_xlim()[1], sub.get_ylim()[1])

    for subax in ax:
        for sub in subax:
            sub.set_xlim(min_xy,max_xy)
            sub.set_ylim(min_xy,max_xy)

    ax[0][0].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[0][1].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[1][0].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)
    ax[1][1].plot((min_xy,max_xy),(min_xy,max_xy),'--',color='0.7',zorder=0)

    ax[0][0].set_xlabel(csv1.split('/')[-1].split('_sorted')[0]+' rep1', fontsize=14)
    ax[0][0].set_ylabel(csv1.split('/')[-1].split('_sorted')[0]+' rep2', fontsize=14)

    ax[0][1].set_xlabel(csv2.split('/')[-1].split('_sorted')[0]+' rep1', fontsize=14)
    ax[0][1].set_ylabel(csv2.split('/')[-1].split('_sorted')[0]+' rep2', fontsize=14)

    ax[1][0].set_xlabel(csv1.split('/')[-1].split('_sorted')[0]+' rep1', fontsize=14)
    ax[1][0].set_ylabel(csv2.split('/')[-1].split('_sorted')[0]+' rep1', fontsize=14)

    ax[1][1].set_xlabel(csv1.split('/')[-1].split('_sorted')[0]+' rep2', fontsize=14)
    ax[1][1].set_ylabel(csv2.split('/')[-1].split('_sorted')[0]+' rep2', fontsize=14)

    fig.tight_layout()
    plt.show()
    plt.clf()
    return fig

def main():
    csv1 = sys.argv[1]
    csv2 = sys.argv[2]
    plot_name = sys.argv[3]
    if len(sys.argv) > 4:
        if sys.argv[4] == '--min_overlap':
            min_overlap = float(sys.argv[5])
    else:
        min_overlap = 0.5
        
    df1 = prep_csv(csv1)
    df2 = prep_csv(csv2)
    
    df1 = compare_peaks_MACS(df1, df2, min_overlap=min_overlap)
    
    fig = scatter_plots(df1, csv1, csv2)
    fig.savefig(plot_name+'.pdf', format='pdf', bbox_inches='tight')
    
if __name__ == "__main__":
    main()
    
    