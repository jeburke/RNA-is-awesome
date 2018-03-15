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
import pysam

def run(cmd, logfile):
    p = subprocess.Popen(cmd, shell=True, universal_newlines=True, stdout=logfile, stderr=logfile)
    ret_code = p.wait()
    logfile.flush()
    return ret_code

def align_fastq_chip(directory, threads=1, organism=None, adaptor='GATCGGAAGA', gff3=None, bowtie_ix=None):
    '''Automatically aligns all fastq.gz files in a drectory using TopHat'''
    if directory[-1] != '/':
        directory = directory+'/'
        
    if organism is None:
        gff3=gff3
        bowtie_ix=bowtie_ix
    elif 'crypto' in organism.lower():
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
    elif 'candida' in organism.lower() or 'albicans' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/C_albicans'
        gff3 = '/home/jordan/GENOMES/C_albicans_SC5314_version_A21-s02-m09-r10_features.gff'
    else:
        print "Organism not recognized"
        return None
        
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
                        for p in processes:
                            p.wait()
                    
                    #ret_code = run(cutadapt_args, logfile)

        for n, fastq in enumerate(fastq_list):
            prefix = fastq.split('/')[-1].split('.fastq')[0]
            print prefix
            trim_name = prefix+'_trim.fastq'
            with open(prefix+'.log', 'a') as logfile:
                if len([x for x in os.listdir(directory+'MULTI/') if prefix in x]) == 0:
                    logfile.write('\n**********Bowtie multi alignment output**********\n')
                    print "Aligning allowing multiple alignments (random assignment)..."
                    
                    multi_args = "bowtie -p{0} -v2 -M1 --best --un {1}_multi_un.fastq --max {1}_multi.fastq {2} -q {3} --sam {1}_multi.sam".format(str(threads), directory+'MULTI/'+prefix, bowtie_ix, trim_name)
                    ret_code = run(multi_args, logfile)
                
                if len([x for x in os.listdir(directory+'UNIQUE/') if prefix in x]) == 0:
                    logfile.write('\n**********Bowtie unique alignment output**********\n')
                    print "Aligning allowing only unique alignments...\n"
                    
                    unique_args = "bowtie -p{0} -v2 -m1 --un {1}_unique_un.fastq  {2} -q {3} --sam {1}_unique.sam".format(str(threads), directory+'UNIQUE/'+prefix, bowtie_ix, trim_name)
                    ret_code = run(unique_args, logfile)
    
    for subdir in ('MULTI/','UNIQUE/'):
        sam_files = [x for x in os.listdir(directory+subdir) if x.endswith('.sam')]
        processes = []
        for n, sam in enumerate(sam_files):
            name = sam.split('.sam')[0]
            args = "samtools view -Sbo {0}.bam {0}.sam".format(directory+subdir+name)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            
            if (n+1) % threads == 0 or n+1 == len(sam_files):
                for p in processes:
                    p.wait()
                    
        bam_files = [x for x in os.listdir(directory+subdir) if x.endswith('bam') and not x.endswith('_sorted.bam')]
        processes = []
        for n, bam in enumerate(bam_files):
            name = bam.split('.bam')[0]
            args = "samtools sort {0}.bam -o {0}_sorted.bam".format(directory+subdir+name)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            if (n+1) % threads == 0 or n+1 == len(bam_files):
                for p in processes:
                    p.wait()
    
        sorted_bams = [x for x in os.listdir(directory+subdir) if x.endswith('_sorted.bam')]
        processes = []
        for n, bam in enumerate(sorted_bams):
            print bam
            args = "samtools index {0}".format(directory+subdir+bam)
            processes.append(subprocess.Popen(args.split(' '), stdout=subprocess.PIPE))
            if (n+1) % threads == 1 or n+1 == len(sorted_bams):
                for p in processes:
                    p.wait()

        # Remove processing intermediates
        for sam in [x for x in os.listdir(directory+subdir) if x.endswith('.sam')]:
            os.remove(directory+subdir+sam)
        for bam in [x for x in os.listdir(directory+subdir) if x.endswith('.bam')]:
            if not bam.endswith('_sorted.bam'):
                os.remove(directory+subdir+bam)
    
    for trim in [x for x in os.listdir(directory) if x.endswith('trim.fastq')]:
        os.remove(directory+trim)

def main():
    alt_adaptor = None
    gff3 = None
    index = None
    for n, arg in enumerate(sys.argv):
        if arg == '-h' or arg == '--help':
            print "\nUsage:\npython ChIP_tools.py --directory fastq_directory --threads num_threads --organism crypto/pombe/cerevisiae/candida <--adaptor GATCGGAAGA> <--gff3 gff3_file> <--index bowtie_index_prefix>\n"
            print "Note: --adaptor argument is optional and will default to the one shown above\n
            --gff3 and --index must both be provided if not using the default files for your organism.\n"
            return None
        
        elif arg == '--directory':
            directory = sys.argv[n+1]
        elif arg == '--threads':
            threads = int(sys.argv[n+1])
        elif arg == '--organism':
            organism = sys.argv[n+1]
        elif arg == '--adaptor':
            alt_adaptor = sys.argv[n+1]
        elif arg == '--gff3':
            gff3 = sys.argv[n+1]
        elif arg == '--index':
            index = sys.argv[n+1]
    
    if alt_adaptor is not None:
        adaptor = alt_adaptor
    else:
        adaptor = 'GATCGGAAGA'
        
    if gff3 is not None and index is not None:
        organism=None
        
    align_fastq_chip(directory, threads=threads, organism=organism, adaptor=adaptor, gff3=gff3, bowtie_ix=index)
    
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
    transcripts = set()
    transcripts_for_df = []
    for ix, r in df.iterrows():
        r_tx = ''
        chrom_tx = {k:v for k, v in tx_dict.items() if v[3] == r['chr']}
        for tx, info in chrom_tx.iteritems():
            if tx[:-2] not in r_tx:
                if info[2] == '+':
                    if r['abs_summit'] in range(info[0]-500,info[1]):
                        r_tx = r_tx+tx[:-2]+','
                        transcripts.add(tx[:-2])
                if info[2] == '-':
                    if r['abs_summit'] in range(info[0],info[1]+500):
                        r_tx = r_tx+tx[:-2]+','
                        transcripts.add(tx[:-2])
        if r_tx == '':
            for tx, info in chrom_tx.iteritems():
                if tx[:-2] not in r_tx:
                    if info[2] == '+':
                        if r['abs_summit'] in range(info[0]-1000,info[1]):
                            r_tx = r_tx+tx[:-2]+','
                            transcripts.add(tx[:-2])
                    if info[2] == '-':
                        if r['abs_summit'] in range(info[0],info[1]+1000):
                            r_tx = r_tx+tx[:-2]+','
                            transcripts.add(tx[:-2])
        transcripts_for_df.append(r_tx)
    df.loc[:,'transcript'] = transcripts_for_df
    return df, transcripts
        
def compare_MACS_output(rep1_xls, rep2_xls, untagged_xls, organism, return_df=False):
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
    enrich2 = []
    starts = []
    ends = []
    for ix, r in df1.iterrows():
        r_rep = False
        peak = range(r['start'],r['end'])
        rep2_matches = df2[(df2['chr'] == r['chr']) & (df2['abs_summit'].isin(peak))]
        if len(rep2_matches) > 0:
            rep.append(True)
            r_rep = True
            
            # Determine the nearest absolute summit and record fold enrichment in replicate
            s = (rep2_matches['abs_summit']-r['abs_summit']).apply(abs)
            min_dist = s[s == min(s)].index[0]
            enrich2.append(rep2_matches.loc[min_dist,'fold_enrichment'])
            
            # Contract peak window to the reproducible region
            new_start = max(rep2_matches.loc[min_dist,'start'], r['start'])
            new_end = min(rep2_matches.loc[min_dist,'end'], r['end'])
            starts.append(new_start)
            ends.append(new_end)
        else:
            rep.append(False)
            enrich2.append(np.NaN)
            starts.append(r['start'])
            ends.append(r['end'])
            
    df1.loc[:,'In replicate'] = rep
    df1.loc[:,'fold_enrichment2'] = enrich2
    df1.loc[:,'start'] = starts
    df1.loc[:,'end'] = ends
    
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
    print "\nPeaks in "+rep1_xls+":"
    print len(df1)
    print "Peaks in "+rep2_xls+":"
    print len(df2)
    filtered = df1[df1['In replicate'] == True]
    print "Peaks in both replicates:"
    print len(filtered)
    filtered = filtered[filtered['In untagged'] == False]
    print "Peaks not in untagged:"
    print len(filtered)
    filtered = filtered.drop(['In replicate','In untagged'], axis=1)
    
    # Remove redundant peaks
    filtered = filtered.sort_values(['fold_enrichment'], ascending=False)
    filtered = filtered.drop_duplicates(subset=['chr','start','end'])
    print "Peaks remaining after removing summits within the same peak region:"
    print len(filtered)
    
    # Add transcripts
    filtered, transcripts = add_transcript(filtered, gff3, organism=organism)
    
    filtered.to_csv(rep1_xls.split('.xls')[0]+'_comparison.csv', index=False)
    
    with open(rep1_xls.split('.xls')[0]+'_genes_with_peaks.txt', 'w') as fout:
        for transcript in transcripts:
            fout.write(transcript+'\n')
    if return_df is True:
        return filtered
    else:
        return None

def wt_v_mut_MACS(df1, df2, min_overlap=0.5):
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
    
    print len(df1)
    counter = 0
    for ix, r in df2.iterrows():
        match = False
        chrom_df = df1[df1['chr'] == r['chr']]
        range1 = set(range(r['start'],r['end']))
        for ix2 in chrom_df.index:
            range2 = range(chrom_df.loc[ix2,'start'], chrom_df.loc[ix2,'end'])
            if len(range1.intersection(range2))/float(len(range1)) >= min_overlap:
                match = True
                break
        if match is False:
            new_row = df2.loc[ix,:]
            new_row = new_row.rename({'fold_enrichment': 'fold_enrichment mut1', 'fold_enrichment2': 'fold_enrichment mut2'})
            df1 = df1.append(new_row)
            df1 = df1.reset_index(drop=True)
            counter += 1
    print counter
    print len(df1)
    
    df1.index = df1['chr'].str.cat(df1['start'].apply(str),sep=':').str.cat(df1['end'].apply(str),sep='-')
            
    return df1    
    
def MACS_scatter_plots(df, xls_pair1, xls_pair2):
    WT_a = [x for x in df.columns if x in xls_pair1[0]][0]
    WT_b = [x for x in df.columns if x in xls_pair1[1]][0]
    MUT_a = [x for x in df.columns if x in xls_pair2[0]][0]
    MUT_b = [x for x in df.columns if x in xls_pair2[1]][0]
    
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10,10))
    ax[0][0].scatter(df[WT_a].apply(np.log2), df[WT_b].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[0][1].scatter(df[MUT_a].apply(np.log2), df[MUT_b].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[1][0].scatter(df[WT_a].apply(np.log2), df[MUT_a].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[1][1].scatter(df[WT_b].apply(np.log2), df[MUT_b].apply(np.log2), s=15, color='0.3', alpha=0.5)

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

    ax[0][0].set_xlabel(xls_pair1[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][0].set_ylabel(xls_pair1[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    ax[0][1].set_xlabel(xls_pair2[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][1].set_ylabel(xls_pair2[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    ax[1][0].set_xlabel(xls_pair1[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[1][0].set_ylabel(xls_pair2[0].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    ax[1][1].set_xlabel(xls_pair1[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)
    ax[1][1].set_ylabel(xls_pair2[1].split('/')[-1].split('_sorted')[0]+'\n Enrichment over WCE', fontsize=14)

    fig.tight_layout()
    plt.show()
    plt.clf()
    return fig

def genome_tiles(chromosome_sizes, tile_size=1000):
    names = []
    chroms = []
    start = []
    end = []
    with open(chromosome_sizes) as f:
        for line in f:
            chrom = line.split('\t')[0]
            size = int(line.split('\t')[1].strip())
            n = 0
            while n*tile_size < size-tile_size:
                chroms.append(chrom)
                names.append(chrom+'-'+str(n))
                start.append(n*tile_size)
                end.append(n*tile_size+tile_size)
                n += 1
            chroms.append(chrom)
            names.append(chrom+'-'+str(n))
            start.append(n*tile_size)
            end.append(size)
    df = pd.DataFrame(index=names)
    df.loc[:,'chromosome'] = chroms
    df.loc[:,'start'] = start
    df.loc[:,'end'] = end
    return df

def count_reads_in_tiles(bam_list, WCE, chromosome_sizes, tile_size=1000):
    df = genome_tiles(chromosome_sizes)
    
    bam_list.append(WCE)
    names = []
    for bam in bam_list:
        bam_name = bam.split('/')[-1].split('_sorted.bam')[0]
        if bam == WCE:
            WCE_name = bam_name
        else:
            names.append(bam_name)
        open_bam = pysam.Samfile(bam)

        read_dens = []
        for ix, r in df.iterrows():
            read_counter = 0
            try:
                for read in open_bam.fetch(r['chromosome'],r['start'],r['end']):
                    read_counter += 1
                read_dens.append(read_counter/((r['end']-r['start'])/1000.))
            except ValueError:
                read_dens.append(np.NaN)
        
        df.loc[:,bam_name] = read_dens
        df.loc[:,bam_name] = df[bam_name].divide(sum(df[bam_name].dropna())/1000000.)
    for name in names:
        df.loc[:,name+' norm to WCE'] = df[name]/df[WCE_name]
    
    df = df.dropna(how='any')
    return df


def MACS_peak_RPKM_scatters(xls_pair1, xls_pair2, untagged_xls, bam_list, WCE_bam_list, name, organism='crypto', min_overlap=0.5, exclude_chrom=None):
    if len(bam_list) != len(WCE_bam_list):
        if len(WCE_bam_list) != 1:
            print "Must provide only 1 WCE bam file or a matched WCE bam file for each sample in the bam_list!"
            return None
    
    # First run the comparison script to get peaks that align between replicates and to contract peak regions
    df_wt = compare_MACS_output(xls_pair1[0], xls_pair1[1], untagged_xls, organism, return_df=True)
    df_mut = compare_MACS_output(xls_pair2[0], xls_pair2[1], untagged_xls, organism, return_df=True)
    
    # Now need to run through and compare peaks between genotypes - see code for MACS enrichment scatters
    compare_df_A = wt_v_mut_MACS(df_wt, df_mut, min_overlap=min_overlap)
    compare_df_A.to_csv(name+'_MACS_wt_mut_merge.csv')
    
    # Generate dictionary of peaks to quantitate
    peak_dict = {}
    for ix, r in compare_df_A.iterrows():
        key = r['chr']+':'+str(r['start'])+'-'+str(r['end'])
        peak_dict[key] = [r['start'],r['end'],'+',r['chr'], r['transcript']]
    
    if exclude_chrom is not None:
        if type(exclude_chrom) == list:
            for chrom in exclude_chrom:
                peak_dict = {k:v for k,v in peak_dict.items() if v[3] != chrom}
        elif type(exclude_chrom) == str:
            print exclude_chrom
            print len(peak_dict)
            peak_dict = {k:v for k,v in peak_dict.items() if v[3] != exclude_chrom}
            print len(peak_dict)
    
    # Then need to count reads in each bam file
    print "\nCalculating RPKM in peak regions..." 
    data_dict = {}
    bam_names = []
    for bam_file in bam_list: 
        bam_name = bam_file.split('_sorted')[0].split('/')[-1]
        bam_names.append(bam_name)
        print bam_name
        data_dict[bam_name] = Compare_RPKM.count_reads_in_ChIP(peak_dict, bam_file)
        new_ix = data_dict[bam_name].index

    # Assemble dataframe from results
    data_df_columns = data_dict.keys()
    data_df = pd.DataFrame(index=new_ix, columns=data_df_columns)
    for col in data_df.columns:
        data_df.loc[:,col] = data_dict[col]
    in_wt = []
    in_mut = []
    for ix, r in data_df.iterrows():
        if str(compare_df_A.loc[ix,'fold_enrichment']) == 'nan':
            in_wt.append(False)
        else:
            in_wt.append(True)
        
        if str(compare_df_A.loc[ix,'fold_enrichment mut1']) == 'nan':
            in_mut.append(False)
        else:
            in_mut.append(True)
    data_df.loc[:,'In WT'] = in_wt
    data_df.loc[:,'In Mut'] = in_mut
        
    # Count whole cell extract and normalize
    for n, bam_file in enumerate(WCE_bam_list):
        WCE = Compare_RPKM.count_reads_in_ChIP(peak_dict, bam_file)
        if len(WCE_bam_list) == 1:
            for bam in bam_names:
                data_df.loc[:,bam] = data_df[bam]/WCE
        elif len(WCE_bam_list) == len(bam_list):
            data_df.loc[:,bam_names[n]] = data_df[bam_names[n]]/WCE

    # Now get the background level and apply adjustment
    print "\n Adjusting for biases..."
    chromosome_sizes = '/home/jordan/GENOMES/crypto_for_bedgraph.genome'
    tile_df = count_reads_in_tiles(bam_list, WCE_bam_list[0], chromosome_sizes)
    
    for bam_name in bam_names:
        if bam_name in xls_pair1[0] or bam_name in xls_pair1[1]:
            peak_min = np.percentile(data_df[data_df['In WT'] == True][bam_name], 0.5)
        elif bam_name in xls_pair2[0] or bam_name in xls_pair2[1]:
            peak_min = np.percentile(data_df[data_df['In Mut'] == True][bam_name], 0.5)
        
        
        s = tile_df.sort_values(bam_name+' norm to WCE').reset_index()[bam_name+' norm to WCE']
        top = max(s[(s <= peak_min)].index)
        slope = stats.linregress(tile_df.sort_values(bam_name).reset_index().iloc[1000:top].index,
            tile_df.sort_values(bam_name).reset_index().iloc[1000:top][bam_name])[0]
        print bam_name
        print slope
        
        data_df.loc[:,bam_name] = data_df[bam_name].divide(slope*1000)
    
    # Add transcripts to spreadsheet
    tx_list = []
    for ix, r in data_df.iterrows():
        tx_list.append(peak_dict[ix][4])
        
    data_df.loc[:,'transcript'] = tx_list
    data_df.to_csv(name+'.csv')
    
    # Make scatter plots - 0v1, 2v3, 0v2, 1v3
    fig = MACS_scatter_plots(data_df, xls_pair1, xls_pair2)
    fig.savefig(name+'.pdf',format='pdf',bbox_inches='tight')
    
    #return data_df 
    
def get_transcripts_from_peak_csv(csv):
    df = pd.read_csv(csv, index_col=0)

    transcripts = []
    for ix, r in df.iterrows():
        if type(r['transcript']) == str:
            r_tx = r['transcript'].split(',')
            transcripts = transcripts + r_tx
    transcripts = set([x for x in transcripts if len(x) > 0])
    print len(transcripts)

    with open(csv.split('.csv')[0]+'_transcripts.txt','w') as fout:
        for tx in transcripts:
            fout.write(tx+'\n')
            
def MACS_Z_score(csv, wt, mut):
    data_df = pd.read_csv(csv, index_col=0)
    wt_names = [x for x in data_df.columns if wt in x]
    mut_names = [x for x in data_df.columns if mut in x]

    s1 = (data_df[mut_names[0]]/data_df[wt_names[0]]).apply(np.log2)
    s1_Z = pd.Series(stats.mstats.zscore(s1), index=data_df.index)

    s2 = (data_df[mut_names[1]]/data_df[wt_names[1]]).apply(np.log2)
    s2_Z = pd.Series(stats.mstats.zscore(s2), index=data_df.index)

    up = set(s1_Z[s1_Z >= 1.8].index).intersection(s2_Z[s2_Z >= 1.8].index)
    print "Genes up in mutant: "+str(len(up))
    down = set(s1_Z[s1_Z <= -1.8].index).intersection(s2_Z[s2_Z <= -1.8].index)
    print "Genes down in mutant: "+str(len(down))
    
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10,10))
    ax[0][0].scatter(data_df[wt_names[0]].apply(np.log2), data_df[wt_names[1]].apply(np.log2), s=15, color='0.3', alpha=0.5)
    ax[0][1].scatter(data_df[mut_names[0]].apply(np.log2), data_df[mut_names[1]].apply(np.log2), s=15, color='0.3', alpha=0.5)

    up_df = data_df[data_df.index.isin(up)]
    up_df.to_csv(csv.split('.csv')[0]+'_up_in_mut.csv')
    down_df = data_df[data_df.index.isin(down)]
    down_df.to_csv(csv.split('.csv')[0]+'_down_in_mut.csv')
    other_df = data_df[(~data_df.index.isin(up)) & (~data_df.index.isin(down))]

    ax[1][0].scatter(other_df[wt_names[0]].apply(np.log2), other_df[mut_names[0]].apply(np.log2), s=15, color='0.7', alpha=0.5)
    ax[1][0].scatter(up_df[wt_names[0]].apply(np.log2), up_df[mut_names[0]].apply(np.log2), s=15, color='orangered')
    ax[1][0].scatter(down_df[wt_names[0]].apply(np.log2), down_df[mut_names[0]].apply(np.log2), s=15, color='royalblue')

    ax[1][1].scatter(other_df[wt_names[1]].apply(np.log2), other_df[mut_names[1]].apply(np.log2), s=15, color='0.7', alpha=0.5)
    ax[1][1].scatter(up_df[wt_names[1]].apply(np.log2), up_df[mut_names[1]].apply(np.log2), s=15, color='orangered')
    ax[1][1].scatter(down_df[wt_names[1]].apply(np.log2), down_df[mut_names[1]].apply(np.log2), s=15, color='royalblue')

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

    ax[0][0].set_xlabel(wt_names[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][0].set_ylabel(wt_names[1]+'\n Enrichment over WCE', fontsize=14)

    ax[0][1].set_xlabel(mut_names[0]+'\n Enrichment over WCE', fontsize=14)
    ax[0][1].set_ylabel(mut_names[1]+'\n Enrichment over WCE', fontsize=14)

    ax[1][0].set_xlabel(wt_names[0]+'\n Enrichment over WCE', fontsize=14)
    ax[1][0].set_ylabel(mut_names[0]+'\n Enrichment over WCE', fontsize=14)

    ax[1][1].set_xlabel(wt_names[1]+'\n Enrichment over WCE', fontsize=14)
    ax[1][1].set_ylabel(mut_names[1]+'\n Enrichment over WCE', fontsize=14)

    fig.tight_layout()
    fig.savefig(csv.split('.csv')[0]+'_Z_scatter.pdf', format='pdf', bbox_inches='tight')
    plt.show()
    plt.clf()
    
def RNAseq_ChIP_overlap_volcano(DESeq2_csv, ChIP_gene_list):
    RNA = GT.load_DESeq2_results([DESeq2_csv])
    chip = pd.read_csv(ChIP_gene_list, index_col=0, header=None)
    RNA = RNA[RNA.index.isin(chip.index)]
    GT.volcano_plot(RNA)
