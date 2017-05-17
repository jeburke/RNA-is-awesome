import sys
import os
import subprocess
import pandas as pd
from sklearn import cluster
from matplotlib import pyplot as plt
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

def align_fastq(directory, threads=1, organism='crypto'):
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
        print column
        new_col_tuples.append((new_name,column))
    new_cols = pd.MultiIndex.from_tuples(new_col_tuples)
    df.columns = new_cols
    return df
    
def load_DESeq2_results(csv_list):
    if csv_list == 'all':
        csv_list = []
        cwd = os.getcwd()
        for file in os.listdir(cwd):
            if file.endswith('.csv'):
                csv_list.append(file)
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
    column_names = []
    if sample_names is not None:
        for name in sample_names:
            column_names.append((name,'log2FoldChange'))
    else:
        for column in dataframe.columns:
            if column[1] == 'log2FoldChange':
                column_names.append(column)
        
    data = dataframe[column_names]
    data = data.dropna(how='any')

    # Convert DataFrame to matrix
    mat = data.as_matrix()

    # Using sklearn
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
    plt.xticks([x+0.5 for x in range(len(data.columns))], column_names, rotation="vertical")
    for size in cluster_sizes:
        ax.plot([0,len(data.columns)-1], [size,size], '-', color='0.5')
        
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(pcm, cax=cax)
    plt.show()
    
    return data

def list_of_genes_in_cluster(data, cluster_index, name=None):
    in_cluster = data[data['cluster'] == cluster_index]
    gene_list = in_cluster.index.tolist()
    print len(gene_list)
    
    if name is None:
        name = 'cluster'+str(cluster_index)
    
    with open(name+'_genes.txt', 'w') as fout:
        for gene in gene_list:
            fout.write(gene+'\n')
    
    return gene_list
        