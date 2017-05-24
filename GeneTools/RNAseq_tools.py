import sys
import os
import subprocess
import pandas as pd
from sklearn import cluster
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib_venn import venn2, venn2_circles
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
import SPTools as SP
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/GeneTools/')
import Annotation_tools
from math import factorial
from decimal import Decimal

font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12,
        }

def align_fastq(directory, threads=1, organism='crypto'):
    '''Automatically aligns all fastq.gz files in a drectory using TopHat'''
    if directory[-1] != '/':
        directory = directory+'/'
        
    if 'crypto' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/Crypto_for_gobs'
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/S288C/S288C'
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        if 'DBK' in organism:
            bowtie_ix = '/home/jordan/GENOMES/S288C/S288C_DBK'
    elif 'pombe' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/POMBE/Spombe'
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        
    fastq_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("fastq.gz"):
            fastq_list.append(directory+file)
    
    ##Check to make sure tophat hasn't already completed on these files.
    for name in os.listdir(directory):
        if os.path.isdir(directory+name):
            for fastq in fastq_list:
                if directory+name == fastq.split('.fastq.gz')[0]:
                    if 'accepted_hits.bam' in os.listdir(directory+name):
                        fastq_list.remove(fastq)
                        print "Tophat alread completed on "+fastq
    
    for fastq in fastq_list:
        prefix = fastq.split('.fastq')[0]
        args = 'tophat2 --read-mismatches 3 --read-gap-length 2 --read-edit-dist 3 --min-anchor-length 8 --splice-mismatches 1 --min-intron-length 20 --max-intron-length 2000 --max-insertion-length 3 --max-deletion-length 3 --num-threads {0} --max-multihits 100 --library-type fr-firststrand --segment-mismatches 3 --no-coverage-search --segment-length 20 --min-coverage-intron 10 --max-coverage-intron 100000 --min-segment-intron 50 --max-segment-intron 500000 --b2-sensitive --bowtie1 -G {1} -o {2} {3} {4}'.format(str(threads), gff3, prefix, bowtie_ix, fastq)
        #print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        p.wait()
        
        print '\n'
        print prefix
	if 'align_summary.txt' in os.listdir(prefix):
            with open(prefix+'/align_summary.txt','r') as f:
                for line in f:
                    print line
                print '\n'
            
    print 'Finished'

def sort_index_bam_files(base_dir):
    '''Finds all TopHat output files and sorts and indexes BAM files'''
    try:
        check_samtools = subprocess.check_output(['samtools','--version'])
        version = check_samtools.split('\n')[0].split(' ')[1]
    except subprocess.CalledProcessError:
        version = 'old'
    print version
    
    
    if base_dir[-1] != '/':
        base_dir = base_dir+'/'
        
    tophat_out = []
    names = []
    for name in os.listdir(base_dir):
        if os.path.isdir(base_dir+name):
            if 'accepted_hits.bam' in os.listdir(base_dir+name):
                if name.split('/')[-1]+'_sorted.bam' not in os.listdir(base_dir):
                    tophat_out.append(base_dir+name+'/accepted_hits.bam')
                    names.append(base_dir+name)
    
    n=0
    for n in range(len(tophat_out)):
        name = names[n]
        args = 'samtools sort {0} -o {1}'.format(tophat_out[n], name+'_sorted.bam')
        if version == 'old':
            args = 'samtools sort {0} {1}'.format(tophat_out[n], name+'_sorted')
        print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        p.wait()
    
    sorted_bams = []
    for file in os.listdir(base_dir):
        if file.endswith('_sorted.bam'):
            if file+'.bai' not in os.listdir(base_dir):
                sorted_bams.append(base_dir+file)
              
    for bam in sorted_bams:
        args = 'samtools index {0}'.format(bam)
        print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        p.wait()
        
def count_with_HTseq(base_dir, organism='crypto'):
    '''Counts reads in every transcript with HTseq'''
    if base_dir[-1] != '/':
        base_dir = base_dir+'/'
    
    if 'crypto' in organism.lower():
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        if 'anti' in organism.lower():
            gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113_anti.gff3'
    elif 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
    
    
    bam_files = []
    for file in os.listdir(base_dir):
        if file.endswith('_sorted.bam'):
            if file.split('_sorted.bam')[0]+'.htseq' not in os.listdir(base_dir):
                bam_files.append(file)
    
    for bam in bam_files:
        args = 'htseq-count -f bam -s reverse -t gene -i ID {0} {1}'.format(base_dir+bam, gff3)
        print args
        
        out = subprocess.check_output(args.split(' '))
        with open(base_dir+bam.split('_sorted')[0]+'.htseq','w') as fout:
            fout.write(out)
        
def main():
    '''Runs the first part of this module to align, sort, index and count reads'''
    # Input will be RNAseq_tools.py directory number_of_threads_tophat organism
    # Current organism options: 'crypt', 'pombe', 'cerevisiae'
    if len(sys.argv) == 3:
        organism = 'crypto'
    else:
        organism = sys.argv[3]
        
    directory = sys.argv[1]
    threads = sys.argv[2]
    
    print "********Aligning reads with Tophat********\n"
    align_fastq(directory, threads=threads, organism=organism)
    print "********Sorting and indexing BAM files********\n"
    sort_index_bam_files(directory)
    print "********Counting reads in transcripts with HTseq********\n"
    count_with_HTseq(directory, organism=organism)
    
if __name__ == "__main__":
    main()

    
def add_col_level(df, new_name):
    new_col_tuples = []
    for column in df:
        #print column
        new_col_tuples.append((new_name,column))
    new_cols = pd.MultiIndex.from_tuples(new_col_tuples)
    df.columns = new_cols
    return df
    
def load_DESeq2_results(csv_list):
    '''This function loads DESeq2 output csv files into a dataframe.
    csv_list can be either a list of the csv files - ["file1.csv","file2.csv"] or "all".'''
    if csv_list == 'all':
        csv_list = []
        cwd = os.getcwd()
        for file in os.listdir(cwd):
            if file.endswith('.csv'):
                csv_list.append(file)
        print csv_list
    
    n=0
    for n, file in enumerate(csv_list):
        if n == 0:
            df = pd.read_csv(file, index_col=0)
            df = add_col_level(df, file.split('/')[-1].split('.csv')[0])
        else:
            new_df = pd.read_csv(file, index_col=0)
            new_df = add_col_level(new_df, file.split('/')[-1].split('.csv')[0])
            df = df.merge(new_df, left_index=True, right_index=True)
    return df

def RNAseq_clustered_heatmap(dataframe, sample_names=None, n_clusters=10):
    '''Generates a clustered heatmap using kmeans clustering and returns a dataframe with cluster indeces.
    If no sample names are provided, clusters across all samples in the dataframe.
    n_clusters is the number of clusters'''
    
    column_names = []
    if sample_names is not None:
        for name in sample_names:
            column_names.append((name,'log2FoldChange'))
    else:
        sample_names = []
        for column in dataframe.columns:
            if column[1] == 'log2FoldChange':
                column_names.append(column)
                sample_names.append(column[0])
        
    data = dataframe[column_names]
    data = data.dropna(how='any')

    # Convert DataFrame to matrix
    mat = data.as_matrix()

    # Using sklearn, cluster data
    km = cluster.KMeans(n_clusters=n_clusters)
    km.fit(mat)

    # Get cluster assignment labels
    labels = km.labels_

    # Format results as a DataFrame
    data['cluster'] = labels
    data = data.sort_values('cluster')
    cluster_sizes = []
    base=0
    for index in set(data['cluster'].tolist()):
        cluster_sizes.append(len(data[data['cluster'] == index])+base)
        base += len(data[data['cluster'] == index])
        
        
    # Get sizes of clusters
    all_data = []
    for column in data[data.columns[:-1]].columns:
        all_data = all_data + data[column].tolist()

    # Plot data as heatmap
    data_min = min(all_data)
    data_max = max(all_data)
    both_max = max(data_min*-1, data_max)
    fig = plt.figure(figsize=(3,12))
    ax = fig.add_subplot(111)
    pcm = ax.pcolor(range(len(data.columns)), range(len(data.index)), data[data.columns[:-1]], cmap='PuOr_r', norm=colors.Normalize(vmin=both_max*-1, vmax=both_max))
    ax.yaxis.set_ticks(cluster_sizes)
    ax.yaxis.set_ticklabels(range(len(cluster_sizes)))
    plt.xticks([x+0.5 for x in range(len(data.columns))], sample_names, rotation="vertical")
    for size in cluster_sizes:
        ax.plot([0,len(data.columns)-1], [size,size], '-', color='0.5')
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(pcm, cax=cax)
    plt.show()
    
    new_df = dataframe.merge(pd.DataFrame(data[('cluster',)]), right_index=True, left_index=True)
    
    return new_df

def list_of_genes_in_cluster(data, cluster_index, name=None, annotate=False, organism=None):
    '''Pull out a cluster from RNAseq_clustered_heatmap.
    Generates a csv file with annotation of those genes.
    data: output from RNAseq_clustered_heatmap
    cluster_index: the number at the top of the cluster on the heatmap display
    name: prefix for the file name
    annotate: True if you want to generate a csv file with annotation
    organism: "pombe" or "crypto" to apply the correct annotation'''
    
    in_cluster = data[data['cluster'] == cluster_index]
    
    if name is None:
        name = 'cluster'+str(cluster_index)
    
    in_cluster.to_csv(name+'_genes.csv')
            
    if annotate is True:
        if organism is None:
            print "Must include the organism!"
        elif 'crypto' in organism.lower():
            Annotation_tools.crypto_annotate(name+'_genes.csv')
        elif 'pombe' in organism.lower():
            Annotation_tools.pombe_annotate(name+'_genes.csv')
        else:
            print "Organism not recognized. Only 'pombe' and 'crypto' supported at this time"

def volcano_plot(df, sample_names=None, annotate=False, organism=None):
    '''Generates volcano plots and csv files of significantly changed transcripts for each sample.
    df: output from load_DESeq2_results
    sample_names: a list of names if you don't want to analyze all the samples
    annotate: True if you want to generate a csv file with annotation
    organism: "pombe" or "crypto" to apply the correct annotation'''
    
    # Find all sample names
    if sample_names is None:
        sample_names = []
        for column in df.columns:
            if column[1] == 'log2FoldChange':
                sample_names.append(column[0])

    for name in sample_names:
        
        ## Find significantly changed transcripts
        sig = df[df[(name,'pvalue')] <= 0.01][name]
        sig = sig[sig['log2FoldChange'] >= 1]
        sig2 = df[df[(name,'pvalue')] <= 0.01][name]
        sig2 = sig2[sig2['log2FoldChange'] <= -1]
        print name
        print "Increased:"
        print len(sig)
        print "Decreased:"
        print len(sig2)
        
        #Build plot
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(df[(name,'log2FoldChange')], df[(name,'pvalue')].apply(np.log).multiply(-1), color='0.5')
        ax.scatter(sig['log2FoldChange'], sig['pvalue'].apply(np.log).multiply(-1), color='orange')
        ax.scatter(sig2['log2FoldChange'], sig2['pvalue'].apply(np.log).multiply(-1), color='darkslateblue')
        xmin = min(df[(name,'log2FoldChange')])-0.1*min(df[(name,'log2FoldChange')])
        xmax = max(df[(name,'log2FoldChange')])-0.1*max(df[(name,'log2FoldChange')])
        x_all = max([xmin*-1,xmax])
        xmin= -1*x_all
        xmax= x_all
        
        ymed = np.percentile(df[(name,'pvalue')].apply(np.log).multiply(-1).dropna(how='any'), 99)
        ymax = ymed + 1.5*ymed
        
        ax.set_ylim([-5,ymax])
        ax.set_xlim([xmin,xmax])
        ax.plot([xmin,xmax],[2,2], '--', color='0.7')
        ax.plot([1,1], [-5,ymax], '--', color='0.7')
        ax.plot([-1,-1], [-5,ymax], '--', color='0.7')
        ax.set_xlabel("log2 fold change")
        ax.set_ylabel("-log10 pvalue")
        ax.set_title(name)
        plt.show()

        fig.savefig(name+'_volcano.pdf',format='pdf')
        plt.clf()
        
        new_sig = sig.append(sig2)
        print len(new_sig)
        new_sig.to_csv(name+'_sig_changed.csv')
        
        if annotate is True:
            if organism is None:
                print "Must include the organism!"
            elif 'crypto' in organism.lower():
                Annotation_tools.crypto_annotate(name+'_sig_changed.csv')
            elif 'pombe' in organism.lower():
                Annotation_tools.pombe_annotate(name+'_sig_changed.csv')
            else:
                print "Organism not recognized. Only 'pombe' and 'crypto' supported at this time"

def NchooseR(n,r):
    ''' Called by gene_venn'''
    c = factorial(n)/Decimal(factorial(r)*factorial(n-r))
    return c

def hypergeometric(N, n, K, J, k):
    '''Called by gene_venn'''
    K = min(K,J)
    p = (NchooseR(K,k)*NchooseR(N-K,n-k))/NchooseR(N,n)
    x = n*k/N
    return p, x

def venn_2sample(n,K,k,J, name1, name2, colors, p, x):
    A_only = K-k
    B_only = J-k
    AB = k
    s = (A_only, B_only, AB)

    v = venn2(subsets=s, set_labels=(name1, name2))

    # Subset colors
    v.get_patch_by_id('10').set_color(colors[0])
    v.get_patch_by_id('01').set_color(colors[1])
    if AB > 0:
        v.get_patch_by_id('11').set_color(colors[2])

    # Subset alphas
    v.get_patch_by_id('10').set_alpha(0.8)
    v.get_patch_by_id('01').set_alpha(0.8)
    if AB > 0:
        v.get_patch_by_id('11').set_alpha(0.9)

    # Border styles
    c = venn2_circles(subsets=s, linestyle='solid')
    for sub_c in c:
        sub_c.set_lw(1.0)       # Line width

    p = '%.1E' % p
    if x >= 5: x = 'Enriched'
    else: x = 'De-enriched'
    
    plt.title(name1+' vs. '+name2+'\n p-value: '+p+', '+x, fontdict=font)
    plt.show()
    plt.savefig(name1+'_'+name2+'_venn.pdf', format='pdf')
    plt.clf()
    
def gene_venn(csv_files, organism):
    '''Finds overlap between 2 or 3 lists of genes.
    csv_files: list of 2 or 3 csv files where the first column is the gene name (make sure the gene name format matches).
    organism: 'crypto', 'cerevisiae' or 'pombe'
    Output: PDF files of venn diagrams (pairwise) and merged csv files containing the overlapping genes.'''
        
    if 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        organism = 'pombe'
    elif 'crypto' in organism.lower() or 'h99' in organism.lower():
        organism = None
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        organism = None
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
    
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    transcripts = tx_dict.keys()
    genes = set([x[:-2] for x in transcripts])
    
    df_dict = {}
    names = []
    for csv in csv_files:
        name = csv.split('/')[-1].split('.')[0]
        df_dict[name] = pd.read_csv(csv, index_col=0)
        df_dict[name] = add_col_level(df_dict[name], name)
        names.append(name)
    
    # N = genome size
    # n = number of genes in analysis (so len(a)+len(b))
    # K = number of genes in group 1 (len(a))
    # k = overlap(len(a&b))
    
    N = len(genes)
    n = len(df_dict[names[0]]) + len(df_dict[names[1]])
    K = len(df_dict[names[0]])
    overlap = set(df_dict[names[0]].index).intersection(df_dict[names[1]].index)
    k = len(overlap)
    J = len(df_dict[names[1]])
    
    p_ab, x_ab = hypergeometric(N,n,K,J,k)
    venn_2sample(n,K,k,J, names[0], names[1], ['crimson','deepskyblue','darkorchid'], p_ab, x_ab)
    
    
    df_ab = df_dict[names[0]].merge(df_dict[names[1]], right_index=True, left_index=True)
    df_ab.to_csv('{0}_{1}_overlap.csv'.format(names[0], names[1]))

    if len(names) == 3:
        n_ac = len(df_dict[names[0]]) + len(df_dict[names[2]])
        overlap_ac = set(df_dict[names[0]].index).intersection(df_dict[names[2]].index)
        k_ac = len(overlap_ac)
        J_ac = len(df_dict[names[2]])
        
        venn_2sample(n_ac,K,k_ac,J_ac, names[0], names[2], ['crimson','gold','darkorange'])
        p_ac = hypergeometric(N,n_ac,K,J_ac,k_ac)
        print "p-value: "+'%.1E' % p_ac
        
        df_ac = df_dict[names[0]].merge(df_dict[names[2]], right_index=True, left_index=True)
        df_ac.to_csv('{0}_{1}_overlap.csv'.format(names[0], names[2]))

        n_bc = len(df_dict[names[1]]) + len(df_dict[names[2]])
        overlap_bc = set(df_dict[names[1]].index).intersection(df_dict[names[2]].index)
        k_bc = len(overlap_bc)
        J_bc = len(df_dict[names[2]])
        K_bc = len(df_dict[names[1]])
        
        venn_2sample(n_bc,K_bc,k_bc,J_bc, names[1], names[2], ['deepskyblue','gold','forestgreen'])
        p_bc = hypergeometric(N,n_bc,K_bc,J_bc,k_bc)
        print "p-value: "+'%.1E' % p_bc
        
        df_bc = df_dict[names[1]].merge(df_dict[names[2]], right_index=True, left_index=True)
        df_bc.to_csv('{0}_{1}_overlap.csv'.format(names[1], names[2]))