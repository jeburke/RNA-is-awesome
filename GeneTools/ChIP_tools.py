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
import subprocess
import os


def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile, stderr=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

def align_fastq_chip(directory, threads=1, organism='crypto', adaptor='GATCGGAAGA'):
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
        if file.lower().endswith("fastq.gz") or file.lower().endswith("fastq"):
            fastq_list.append(directory+file)
    
    if not os.path.exists(directory+'MULTI/'):
        os.makedirs(directory+'MULTI/')
    if not os.path.exists(directory+'UNIQUE/'):
        os.makedirs(directory+'UNIQUE/')
    
    ##Check to make sure tophat hasn't already completed on these files. Do some cleanup
    for fastq in fastq_list:
        bam = None
        sorted_bam = None
        sam = None
        for align_type in ('MULTI/','UNIQUE/'):
            for name in os.listdir(directory+align_type):    
                if fastq.split('.fastq')[0] in name and name.endswith('_sorted.bam'):
                    fastq_list.remove(directory+fastq)
                    print "Bowtie and Samtools view already completed on "+fastq
                    sorted_bam = name
                elif fastq.split('.fastq')[0] in name and name.endswith('.sam'):
                    fastq_list.remove(directory+fastq)
                    print "Bowtie already completed on "+fastq
                    sam = name
                elif fastq.split('.fastq')[0] in name and name.endswith('.bam'):
                    fastq_list.remove(directory+fastq)
                    print "Bowtie already completed on "+fastq
                    bam = name
            if sorted_bam is not None:
                if sam is not None:
                    os.remove(directory+align_type+sam)
                if bam is not None:
                    os.remove(directory+align_type+bam)
            elif sorted_bam is None and bam is not None:
                if sam is not None:
                    os.remove(directory+align_type+sam)
    
    if len(fastq_list) > 0:
        processes = []
        for n, fastq in enumerate(fastq_list):
            prefix = fastq.split('/')[-1].split('.fastq')[0]
            print prefix
            with open(prefix+'.log', 'w') as logfile:
                trim_name = fastq.split('/')[-1].split('.fastq')[0]+'_trim.fastq'
                if trim_name not in os.listdir(directory):
                    print "Trimming reads..."
                    cutadapt_args = "cutadapt -a {0} -m 20 -o {1}_trim.fastq {2}".format(adaptor, prefix, fastq)
                    logfile.write('**********Cutadapt output**********\n')
                    processes.append(subprocess.Popen(cutadapt_args, shell=True, universal_newlines=True, stdout=logfile, stderr=logfile))
                    if (n+1) % threads == 0 or n+1 == len(fastq_list):
                        print "waiting..."
                        for p in processes:
                            p.wait()
                    
                    #ret_code = run(cutadapt_args, logfile)

        for n, fastq in enumerate(fastq_list):
            prefix = fastq.split('/')[-1].split('.fastq')[0]
            print prefix
            with open(prefix+'.log', 'a') as logfile:
                if len([x for x in os.listdir(directory+'MULTI/') if prefix in x]) == 0:
                    print "Aligning allowing multiple alignments (random assignment)..."
                    multi_args = "bowtie -p{0} -v2 -M1 --best --un {1}_multi_un.fastq --max {1}_multi.fastq {2} -q {3} --sam {1}_multi.sam".format(str(threads), directory+'MULTI/'+prefix, bowtie_ix, trim_name)
                    logfile.write('\n**********Bowtie multi alignment output**********\n')
                    ret_code = run(multi_args, logfile)
                
                if len([x for x in os.listdir(directory+'UNIQUE/') if prefix in x]) == 0:
                    print "Aligning allowing only unique alignments...\n"
                    unique_args = "bowtie -p{0} -v2 -m1 --un {1}_unique_un.fastq  {2} -q {3} --sam {1}_unique.sam".format(str(threads), directory+'UNIQUE/'+prefix, bowtie_ix, trim_name)
                    logfile.write('\n**********Bowtie unique alignment output**********\n')
                    ret_code = run(unique_args, logfile)
    
    for subdir in ('MULTI/','UNIQUE/'):
        sam_files = [x for x in os.listdir(directory+subdir) if x.endswith('.sam')]
        print sam_files
        processes = []
        for n, sam in enumerate(sam_files):
            name = sam.split('.sam')[0]
            args = "samtools view -Sbo {0}.bam {0}.sam".format(directory+subdir+name)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            
            if (n+1) % threads == 0 or n+1 == len(sam_files):
                print "waiting..."
                for p in processes:
                    p.wait()
                    
        bam_files = [x for x in os.listdir(directory+subdir) if x.endswith('bam') and not x.endswith('_sorted.bam')]
        print bam_files
        processes = []
        for n, bam in enumerate(bam_files):
            name = bam.split('.bam')[0]
            args = "samtools sort {0}.bam -o {0}_sorted.bam".format(directory+subdir+name)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            if (n+1) % threads == 0 or n+1 == len(bam_files):
                print "waiting..."
                for p in processes:
                    p.wait()
    
        sorted_bams = [x for x in os.listdir(directory+subdir) if x.endswith('_sorted.bam')]
        print sorted_bams
        processes = []
        for n, bam in enumerate(sorted_bams):
            args = "samtools index {0}".format(directory+subdir+bam)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            if (n+1) % threads == 1 or n+1 == len(sorted_bams):
                for p in processes:
                    p.wait()

        # Remove processing intermediates
        for sam in [x for x in os.listdir(directory+subdir) if x.endswith('.sam')]:
            os.remove(directory+subdir+sam)
        for bam in [x for x in os.listdir(directory+subdir) if x.endswith('.bam')]:
            if not name.endswith('_sorted.bam'):
                os.remove(directory+subdir+bam)
    
    for trim in [x for x in os.listdir(directory) if x.endswith('trim.fastq')]:
        os.remove(directory+trim)

def main():
    alt_adaptor = None
    for n, arg in enumerate(sys.argv):
        if arg == '-h' or arg == '--help':
            print "\nUsage:\npython ChIP_tools.py --directory fastq_directory --threads num_threads --organism crypto/pombe/cerevisiae <--adaptor GATCGGAAGA>\n"
            print "Note: --adaptor argument is optional and will default to the one shown above\n"
            return None
        
        elif arg == '--directory':
            directory = sys.argv[n+1]
        elif arg == '--threads':
            threads = int(sys.argv[n+1])
        elif arg == '--organism':
            organism = sys.argv[n+1]
        elif arg == '--adaptor':
            adaptor = sys.argv[n+1]
    
    if alt_adaptor is not None:
        adaptor = alt_adaptor
    else:
        adaptor = 'GATCGGAAGA'
        
    align_fastq_chip(directory, threads=threads, organism=organism, adaptor=adaptor)
    
if __name__ == "__main__":
    main()
        
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
        
def add_transcript(df, gff3, organism=None):
    tx_dict = GT.build_transcript_dict(gff3, organism=organism)
    transcripts = []
    for ix, r in df.iterrows():
        r_tx = ''
        chrom_tx = {k:v for k, v in tx_dict.items() if v[3] == r['chr']}
        for tx, info in chrom_tx.iteritems():
            if tx[:-2] not in r_tx:
                if info[2] == '+':
                    if r['abs_summit'] in range(info[0]-500,info[1]):
                        r_tx = r_tx+tx[:-2]+','
                if info[2] == '-':
                    if r['abs_summit'] in range(info[0],info[1]+500):
                        r_tx = r_tx+tx[:-2]+','
        transcripts.append(r_tx)
    df.loc[:,'transcript'] = transcripts
    df.loc[df[df['transcript'] == ''].index,'transcript'] = None
    return df
        
def compare_MACS_output(rep1_xls, rep2_xls, untagged_xls, organism):
    if 'crypto' in organism.lower():
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
        organism = None
    elif 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        organism = 'pombe'
    elif 'cerev' in organism.lower() or '288' in organism:
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
        organism = None
        
    df1 = pd.read_csv(rep1_xls, sep='\t', skiprows=28, header=0)
    df2 = pd.read_csv(rep2_xls, sep='\t', skiprows=28, header=0)
    df_un = pd.read_csv(untagged_xls, sep='\t', skiprows=28, header=0)
    
    # Determine if peak is in each replicate
    rep = []
    for ix, r in df1.iterrows():
        peak = range(r['start'],r['end'])
        rep2_matches = df2[(df2['chr'] == r['chr']) & (df2['abs_summit'].isin(peak))]
        if len(rep2_matches) > 0:
            rep.append(True)
        else:
            rep.append(False)
    df1.loc[:,'In replicate'] = rep
    
    # Determine if peak is in untagged
    un = []
    for ix, r in df1.iterrows():
        peak = range(r['start'],r['end'])
        un_matches = df_un[(df_un['chr'] == r['chr']) & (df_un['abs_summit'].isin(peak))]
        if len(un_matches) > 0:
            un.append(True)
        else:
            un.append(False)
    df1.loc[:,'In untagged'] = un
    
    # Filter based on reproducibility and untagged
    filtered = df1[(df1['In replicate'] == True) & (df1['In untagged'] == False)]
    filtered = filtered.drop(['In replicate','In untagged'], axis=1)
    
    # Add transcripts
    filtered = add_transcript(filtered, gff3, organism=organism)
    
    filtered.to_csv(rep1_xls.split('.xls')[0]+'_comparison.csv', index=False)