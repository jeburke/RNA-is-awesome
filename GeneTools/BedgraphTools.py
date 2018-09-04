import sys
import os
from subprocess import check_output
from subprocess import call
from subprocess import Popen
from multiprocessing import Pool
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneUtility
import SPTools as SP
from collections import OrderedDict
import csv


#Little function to sum over a window from a series from my sorted-by-gene bedgraph dictionary
def window_count(series, coord_1, coord_2):
    window = series[series.index.isin(range(coord_1,coord_2))]
    total = sum(window)
    return total

def all_windows(stop_codon, end, strand, series, window_size=10):
    # Calculate a window just before the stop codon and then windows of window_size after the stop codon. 
    # Returns list of windows normalized to region upstream of stop codon.
    if strand == '+':
        a = window_count(series, stop_codon-200, stop_codon)/20.
        b = []
        n = stop_codon
        while n < end-window_size:
            b.append(window_count(series, n, n+window_size))
            n += window_size
    if strand == '-':
        a = window_count(series, stop_codon, stop_codon+200)/20.
        b = []
        n = stop_codon
        while n > end+window_size:
            b.append(window_count(series, n-window_size, n))
            n = n-window_size
    if sum(b) < 5000:
        return None
    
    b = [x/float(a) for x in b]
    return b

def pA_windows(stop_codon, end, strand, series, window_size=15):
    if strand == '+':
        b = []
        n = stop_codon
        while n < end-window_size:
            b.append(window_count(series, n, n+window_size))
            n += window_size
    if strand == '-':
        b = []
        n = stop_codon
        while n > end+window_size:
            b.append(window_count(series, n-window_size, n))
            n = n-window_size
    b = [float(x)/(sum(b)*2) for x in b]
    x = [x*(0.1*window_size) for x in range(len(b))]
    return x,b


# Big function for finding genes where (in both replicates) the read density falls off faster in mutant than in wild type.
# Creates plots for each. Can include polyA bedgraphs which will be plotted as bars, but not recommended.
def UTR_windows(tx_dict, wt_reps_plus, wt_reps_minus, mut_reps_plus, mut_reps_minus, window_size=10, wt_pA_bedgraphs=None, mut_pA_bedgraphs=None, cutoff=0.3):
    change_set = set()
    for tx, info in tx_dict.iteritems():
        wt1 = None
        #Calculate normalized windows for all four samples (2 WT replicates and 2 mut replicates)
        strand = info[2]
        if info[2] == '+':
            if tx in wt_reps_minus[0] and tx in wt_reps_minus[1] and tx in mut_reps_minus[0] and tx in mut_reps_minus[1]:
                try:
                    stop_codon = max(info[5])
                except ValueError:
                    stop_codon = info[6]
                #end = info[1]
                end = stop_codon+201
                wt1 = all_windows(stop_codon, end, info[2], wt_reps_minus[0][tx], window_size=window_size)
                wt2 = all_windows(stop_codon, end, info[2], wt_reps_minus[1][tx], window_size=window_size)
                mut1 = all_windows(stop_codon, end, info[2], mut_reps_minus[0][tx], window_size=window_size)
                mut2 = all_windows(stop_codon, end, info[2], mut_reps_minus[1][tx], window_size=window_size)
        elif info[2] == '-':
            if tx in wt_reps_plus[0] and tx in wt_reps_plus[1] and tx in mut_reps_plus[0] and tx in mut_reps_plus[1]:
                try:
                    stop_codon = min(info[4])
                except ValueError:
                    stop_codon = info[6]
                #end = info[0]
                end = stop_codon-200
                wt1 = all_windows(stop_codon, end, info[2], wt_reps_plus[0][tx], window_size=window_size)
                wt2 = all_windows(stop_codon, end, info[2], wt_reps_plus[1][tx], window_size=window_size)
                mut1 = all_windows(stop_codon, end, info[2], mut_reps_plus[0][tx], window_size=window_size)
                mut2 = all_windows(stop_codon, end, info[2], mut_reps_plus[1][tx], window_size=window_size)

        #Check for dropoff from ORF levels in all samples and score
        if wt1 is not None and wt2 is not None and mut1 is not None and mut2 is not None:
            #Add something to check that windows are consecutive
            
            #wt1_score = len([x for x in wt1 if x < 0.3])
            #wt2_score = len([x for x in wt2 if x < 0.3])
            #mut1_score = len([x for x in mut1 if x < 0.3])
            #mut2_score = len([x for x in mut2 if x < 0.3])
            
            wt1_score = 0
            wt2_score = 0
            mut1_score = 0
            mut2_score = 0
            for n in range(len(wt1)-2):
                if wt1[n] < cutoff and wt1[n+1] < cutoff and wt1[n+2] < cutoff: wt1_score += 1
                if wt2[n] < cutoff and wt2[n+1] < cutoff and wt2[n+2] < cutoff: wt2_score += 1
                if mut1[n] < cutoff and mut1[n+1] < cutoff and mut1[n+2] < cutoff: mut1_score += 1
                if mut2[n] < cutoff and mut2[n+1] < cutoff and mut2[n+2] < cutoff: mut2_score += 1
                    
            #If mutant score is higher than wt score, then make a plot of the UTR and add the transcript to a list
            if mut1_score-wt1_score >= 2 and mut2_score-wt2_score >= 2:
                wt1_avg = sum(wt1)/float(len(wt1))
                wt2_avg = sum(wt2)/float(len(wt2))
                mut1_avg = sum(mut1)/float(len(mut1))
                mut2_avg = sum(mut2)/float(len(mut2))
                if mut1_avg*1.1 < wt1_avg and mut2_avg*1.1 < wt2_avg:
                    change_set.add(tx)
                    print tx[:-2]
                    print sum(wt1)/float(len(wt1))
                    print sum(wt2)/float(len(wt2))
                    print sum(mut1)/float(len(mut1))
                    print sum(mut2)/float(len(mut2))

                    fig = plt.figure(figsize=(8, 6), dpi=600)
                    ax = fig.add_subplot(111)
                    
                    ax.plot(range(0, len(wt1)), wt1, color='blue', label='WT rep 1')
                    ax.plot(range(0, len(wt2)), wt2, color='skyblue', label='WT rep 2')
                    ax.plot(range(0, len(mut1)), mut1, color='red', label='Mutant rep 1')
                    ax.plot(range(0, len(mut2)), mut2, color='orange', label='Mutant rep 2')
                    
                    if wt_pA_bedgraphs is not None and mut_pA_bedgraphs is not None:
                        try:
                            if strand == '+':
                                wt_pA = pA_windows(stop_codon, stop_codon+200, strand, wt_pA_bedgraphs[0][tx])
                            if strand == '-':
                                wt_pA = pA_windows(stop_codon, stop_codon-200, strand, wt_pA_bedgraphs[1][tx])
                            #print wt_pA
                            plt.bar(wt_pA[0], wt_pA[1], 1.5, color='blue', alpha=0.6, label='WT polyA sites')
                        except KeyError:
                            print "WT: "+tx
                            
                        try:
                            if strand == '+':
                                mut_pA = pA_windows(stop_codon, stop_codon+200, strand, mut_pA_bedgraphs[0][tx])
                            if strand == '-':
                                mut_pA = pA_windows(stop_codon, stop_codon-200, strand, mut_pA_bedgraphs[1][tx])
                            plt.bar(mut_pA[0], mut_pA[1], 1.5, color='crimson', alpha=0.5, label='Mutant polyA sites')
                        except KeyError:
                            print "Mut: "+tx

                    
                    plt.title(tx[:-2], fontsize=18)
                    plt.xlabel('Distance from stop codon (10 bp increments)', fontsize=16)
                    plt.ylabel('Proportion of reads', fontsize=16)
                    #plt.legend()
                    plt.show()

                    fig.savefig('{0}_falloff.pdf'.format(tx[:-2]), orientation='landscape', format='pdf',transparent=True, frameon=False, bbox_inches='tight', pad_inches=0.5)
                    print "WT1 score = "+str(wt1_score)
                    print "WT2 score = "+str(wt2_score)
                    print "mut1 score = "+str(mut1_score)
                    print "mut2 score = "+str(mut2_score)

    return change_set

def count_aligned_reads(bam):
    command = 'samtools view -F 0x904 -c {0}'.format(bam)
    aligned_reads = check_output(command.split(), shell=False)
    aligned_reads = float(aligned_reads)/1000000
    print '\n'+bam
    print "{0} million reads".format("%.2f" % aligned_reads)
    return {bam:1/aligned_reads}

def run_genomeCoverageBed(command):
    with open('error.log','a') as fout:
        process = Popen(command, shell=True, stdout=fout, stderr=fout)
    ret_code = process.wait()
    return ret_code

def generate_scaled_bedgraphs2(directory, untagged, organism='crypto', start_only=False, stranded=False, threads=1, file_provided=False):
    if 'crypto' in organism.lower():
        genome = '/home/jordan/GENOMES/crypto_for_bedgraph.genome'
    elif 'cerev' in organism.lower():
        genome = '/home/jordan/GENOMES/S288C/S288C_for_bedgraph.genome'
    elif 'pombe' in organism.lower():
        genome = '/home/jordan/GENOMES/POMBE/Sp_for_bg.genome'
    elif 'albicans' in organism.lower() or 'candida' in organism.lower():
        genome = '/home/jordan/GENOMES/C_albicans_for_bg.genome'
    
    bam_list = []
    if not file_provided:
        for file in os.listdir(directory):
            if file.lower().endswith("sorted.bam"):
                bam_list.append(directory+file)
    else:
        bam_list.append(directory)
        base_dir = directory.split('/')[:-1]
        base_dir = '/'.join(base_dir)+'/'
        untagged = [base_dir+x for x in os.listdir(base_dir) if untagged in x and x.endswith('sorted.bam')][0]
        bam_list.append(untagged)
    
    p = Pool(threads)
    totals = {}
    entries = p.map(count_aligned_reads, bam_list)
    for x in entries:
        totals.update(x)
    
    commands = {}
    for n, bam in enumerate(bam_list):
        command = 'genomeCoverageBed -ibam {0} -g {1} -bga -scale {2} '.format(bam, genome, "%.2f" % totals[bam])
        commands[bam] = command
    
    if start_only is True:
        commands = [x+'-5 ' for x in comands]
    if stranded is True:
        stranded_cmds = {}
        for bam, command in commands.iteritems():
            stranded_cmds[bam] = []
            stranded_cmds[bam].append(command+'-strand + ')
            stranded_cmds[bam].append(command+'-strand - ')
        commands = stranded_cmds
    
    final_cmds = []
    for bam, command in commands.iteritems():
        if type(command) == list:
            final_cmds.append(command[0]+'> {0}_plus.bedgraph'.format(bam.split('.bam')[0]))
            final_cmds.append(command[1]+'> {0}_minus.bedgraph'.format(bam.split('.bam')[0]))
        else:
            final_cmds.append(command+'> {0}.bedgraph'.format(bam.split('.bam')[0]))
    
    #for command in final_cmds:
    #    print command
    
    p = Pool(threads)
    codes = p.map(run_genomeCoverageBed, final_cmds)
    
    return codes
        

def generate_scaled_bedgraphs(directory, organism='crypto', start_only=False, stranded=False, file_provided=False):
    if 'crypto' in organism.lower():
        genome = '/home/jordan/GENOMES/crypto_for_bedgraph.genome'
    elif 'cerev' in organism.lower():
        genome = '/home/jordan/GENOMES/S288C/S288C_for_bedgraph.genome'
    elif 'pombe' in organism.lower():
        genome = '/home/jordan/GENOMES/POMBE/Sp_for_bg.genome'
    elif 'albicans' in organism.lower() or 'candida' in organism.lower():
        genome = '/home/jordan/GENOMES/C_albicans_for_bg.genome'
    
    bam_list = []
    if not file_provided:
        for file in os.listdir(directory):
            if file.lower().endswith("sorted.bam"):
                bam_list.append(directory+file)
    else:
        bam_list.append(directory)
            
    total_aligned = []
    for bam in bam_list:
        print bam
        command = 'samtools view -F 0x904 -c {0}'.format(bam)
        aligned_reads = check_output(command.split(), shell=False)
        total_aligned.append(aligned_reads)
        print "Total aligned reads in "+bam
        print aligned_reads
    
    #ratio_list = []
    #n=0
    #for n in range(len(bam_list)):
    #    ratio = float(total_aligned[n])/float(total_aligned[0])
    #    ratio_list.append(1/ratio)
     
    total_aligned = [float(x)/1000000 for x in total_aligned]
    
    for n in range(len(bam_list)):
        out = bam_list[n].split('/')[-1].split('.')[0]
        if stranded is True:
            if start_only is False:
                command1 = 'genomeCoverageBed -ibam {0} -g {1} -bga -strand + -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            elif start_only is True:
                command1 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand + -5 -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            print command1
            bg1 = check_output(command1.split(), shell=False)
            with open('{0}_plus.bedgraph'.format(out),'w') as fout:
                fout.write(bg1)
            if start_only is False:
                command2 = 'genomeCoverageBed -ibam {0} -g {1} -bga -strand - -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            elif start_only is True:
                command2 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand - -5 -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            bg2 = check_output(command2.split(), shell=False)
            with open('{0}_minus.bedgraph'.format(out),'w') as fout:
                fout.write(bg2)
        else:
            if start_only is False:
                command = 'genomeCoverageBed -ibam {0} -g {1} -bga -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            else:
                command = 'genomeCoverageBed -ibam {0} -g {1} -bg -5 -scale {2}'.format(bam_list[n], genome, str(total_aligned[n]))
            print command
            bg = check_output(command.split(), shell=False)
            with open('{0}.bedgraph'.format(out),'w') as fout:
                fout.write(bg)            

def list_bedgraphs(directory):
    plus_list = []
    minus_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("plus.bedgraph"):
            plus_list.append(directory+file)
        elif file.lower().endswith("minus.bedgraph"):
            minus_list.append(directory+file)
    plus_list.sort()
    minus_list.sort()
    bedgraphs = zip(plus_list,minus_list)
    return bedgraphs

############################################################
## Read bedgraph and sort by transcript into dictionary   ##
############################################################

def build_bedgraph_dict(transcript_dict, bedgraph_file):
    '''Function for sorting bedgraph files by gene.
    
    Parameters
    ----------
    bedgraph_file : str
                    Bedgraph file
    transcript_dict : dict
                    Transcript dict generated by build_transcript_dict (in SeqTools module)
                    
    Output
    -------
    sorted bedgraph : file
            Bedgraph file sorted by gene (ends in _by_gene.bedgraph).'''
    
    print datetime.now()
    bedgraph_dict = {}
    transcript_by_chr = {}
    for transcript, coords in transcript_dict.iteritems():
        chromosome = coords[3]
        bedgraph_dict[transcript] = [[],[]]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(transcript)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(transcript)
    
    with open(bedgraph_file, "r") as bedgraph:
        for line in bedgraph:
            columns = re.split(r'\t', line)
            bed_chr = columns[0].strip()
            rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
            if bed_chr in rom_lat:
                bed_chr = rom_lat[bed_chr]
            bed_position = int(columns[1])
            bed_peak = float(columns[3])
            
            if bed_chr in transcript_by_chr:
                transcript_list = transcript_by_chr[bed_chr]
            for transcript in transcript_list:  
                
                #Dictionary for bedgraph. Values will be [list of genomic positions][reads starting at that position]
                if bed_chr == transcript_dict[transcript][3].strip() and bed_position > transcript_dict[transcript][0] and bed_position < transcript_dict[transcript][1]:
                    bedgraph_dict[transcript][0].append(bed_position)
                    bedgraph_dict[transcript][1].append(bed_peak)
   
    with open("{0}_by_gene.bedgraph".format(bedgraph_file.split("/")[-1].split(".")[0]), "a") as fout:
        for transcript, values in bedgraph_dict.iteritems():
            fout.write(transcript+"\n")
            coord_list = map(str, bedgraph_dict[transcript][0])
            coord_line = "\t".join(coord_list)
            fout.write(coord_line+"\n")
            count_list = map(str, bedgraph_dict[transcript][1])
            count_line = "\t".join(count_list)
            fout.write(count_line+"\n")
                                
    bedgraph_dict = collections.OrderedDict(sorted(bedgraph_dict.items()))

    print datetime.now()

def read_sorted_bedgraph(bedgraph_dict_output, transcript_dict, organism=None):
    '''Function for loading sorted bedgraph files.
    
    Parameters
    ----------
    bedgraph_dict_output : str
                    File output by build_bedgraph_dict (ends in _by_gene.bedgraph)
    transcript_dict : dict
                    Transcript dict generated by build_transcript_dict (in SeqTools module)
    organism : str, default ``None``
                    change to 'pombe' if working with S. pombe
                    
    Returns
    -------
    bg_dict : dict
            Dictionary where keys are transcript name and values are pandas series of bedgraph data'''
    
    bg_dict = {}
    count = 0
    dtype = [('coord', int), ('height', float)]
    with open(bedgraph_dict_output,'r') as f:
        n = 0
        for line in f:
            n += 1
            if len(line) > 1:
                
                #Read the transcript line
                if n%3 == 1:
                    tx = line.strip()
                    count += 1
                    #print count
                    if tx[-2] != 'T' and organism != 'pombe':
                        tx = tx+'T0'
                
                #Read the coordinate line
                elif n%3 == 2:
                    coords  = map(int, line.strip().split('\t'))
                
                #Read the values line
                elif n%3 == 0:
                    heights = map(float, line.strip().split('\t'))

                    if tx not in transcript_dict:
                        pass
                    else:
                        all_coords = set(range(min(coords),max(coords)))
                        missing = all_coords.difference(coords)
                        coords = coords + list(missing)

                        #Fill in missing coordinates with zeros
                        zero_fill = [0]*len(missing)
                        heights = heights + zero_fill

                        #Create a pandas series with all coordinates and sort so zeros are inserted appropriately
                        entry = pd.Series(heights, index=coords)
                        entry.sort_index(inplace=True)

                        selected_range = range(transcript_dict[tx][0],transcript_dict[tx][1])
                        entry = entry[entry.index.isin(selected_range)]
                        
                        bg_dict[tx] = entry
    return bg_dict

def decollapse_bedgraph(bedgraph):
    new_bedgraph = bedgraph.split('.bedgraph')[0]+'_full.bedgraph'
    counter1 = 0
    counter2 = 0
    with open(bedgraph) as f:
        with open(new_bedgraph,'w') as fout:
            for line in f:
                counter1 += 1
                data = line.split('\t')
                chrom = data[0]
                start = int(data[1])
                end = int(data[2])
                value = data[3].strip()

                if end-start != 1:
                    n_lines = end-start
                    for n in range(n_lines):
                        counter2 += 1
                        new_start = str(start+n)
                        new_end = str(start+n+1)
                        new_line = '\t'.join([chrom, new_start, new_end, value+'\n'])
                        fout.write(new_line)
                else:
                    counter2 += 1
                    fout.write(line)
    #print counter1
    #print counter2
    
def collapse_bedgraph(bedgraph):
    new_bedgraph = bedgraph+'.tmp'
    counter1 = 0
    counter2 = 0
    with open(bedgraph) as f:
        with open(new_bedgraph, 'w') as fout:
            for line in f:
                data = line.split('\t')
                chrom = data[0]
                start = data[1]
                end = data[2]
                value = data[3].strip()
                if counter1 == 0:
                    prev_chrom = chrom
                    prev_value = value
                    block_start = start
                    prev_end = end
                else:
                    if prev_chrom != chrom:
                        new_line = [prev_chrom, block_start, prev_end, prev_value+'\n']
                        new_line = '\t'.join(new_line)
                        fout.write(new_line)
                        counter2 += 1
                        
                        prev_chrom = chrom
                        prev_value = value
                        prev_end = end
                        block_start = start
                    elif prev_chrom == chrom and prev_value == value:
                        pass
                    elif prev_chrom == chrom and prev_value != value:
                        new_line = [chrom, block_start, end, prev_value+'\n']
                        new_line = '\t'.join(new_line)
                        fout.write(new_line)
                        counter2 += 1
                        
                        prev_value = value
                        block_start = start
                        prev_end = end

                counter1 += 1
    print counter1
    print counter2
    
    os.remove(bedgraph)
    os.rename(bedgraph+'.tmp', bedgraph)
    os.remove(bedgraph+'.tmp')

def bedgraph_reader(bedgraph, chromosomes=None):
    df = pd.read_csv(bedgraph, sep='\t', header=None, names=['chromosome','start','end','RPM'])
    
    if chromosomes is not None:
        df = df[df['chromosome'].isin(chromosomes)]
        
    df.index = df['chromosome'].str.cat(df['start'].apply(str),sep=':')
    
    return df

def write_bedgraph(dataframe, name):
    dataframe.to_csv(name+'combined.bedgraph', index=False, header=False, sep='\t')
##  
def combine_stranded_bedgraph(directory, file_provided=False):
    bg_pairs = []
    if not file_provided:
        for file in os.listdir(directory):
            if file.endswith('plus.bedgraph'):
                bg_pairs.append((file,file.split('plus.bedgraph')[0]+'minus.bedgraph'))
    else:
        name1 = directory.split('.bam')[0]+'plus.bedgraph'
        name2 = directory.split('.bam')[0]+'minus.bedgraph'
        bg_pairs.append((name1, name2))

    for pair in bg_pairs:
        name = pair[0].split('.bedgraph')[0].split('plus')[0]
        if not name.endswith('_'):
            name = name+'_'
        plus = bedgraph_reader(pair[0])
        minus = bedgraph_reader(pair[1])
        minus[3] = minus[3].multiply(-1)

        new = plus.append(minus)
        new = new.sort_values([0,1])
        write_bedgraph(new, name)
##
 
def normalize_bedgraph(tagged, untagged, smooth=False):
    tagged_RPM = bedgraph_reader(tagged)
    untagged_RPM = bedgraph_reader(untagged)
    
    total = tagged_RPM.merge(untagged_RPM, right_index=True, left_index=True, how='left')
    total.loc[:,'norm RPM'] = total['RPM_x']/total['RPM_y']
    
    normalized = total[['chromosome_x','start_x','end_x','norm RPM']]
    normalized = normalized.replace([np.inf,np.inf*-1],np.NaN).dropna(how='any')
    
    normalized.to_csv(tagged.split('.bedgraph')[0]+'_norm.bedgraph', sep='\t', index=False, header=False)
    
    if smooth is False:
        collapse_bedgraph(tagged.split('.bedgraph')[0]+'_norm.bedgraph')
    
def smooth_bedgraphs(bedgraph_list, window):
    for bedgraph in bedgraph_list:
        bg_df = bedgraph_reader(bedgraph)
        new_bg = pd.DataFrame(columns=bg_df.columns)
        for chrom in set(bg_df['chromosome']):
            chrom_df = bg_df[bg_df['chromosome'] == chrom]
            chrom_df = chrom_df.sort_values(['start'])
            new_intensities = chrom_df['RPM'].rolling(window=window, center=True).mean()
            chrom_df.loc[:,'RPM'] = new_intensities
            new_bg = new_bg.append(chrom_df.dropna(how='any'))
            
        new_bg = new_bg.sort_values(['chromosome','start'])
        new_bg.loc[:,'start'] = new_bg['start'].apply(int)
        new_bg.loc[:,'end'] = new_bg['end'].apply(int)
        new_bg.to_csv(bedgraph.split('.bedgraph')[0]+'_{0}bp_smooth.bedgraph'.format(str(window)), sep='\t', index=False, header=False)
        
        collapse_bedgraph(bedgraph.split('.bedgraph')[0]+'_{0}bp_smooth.bedgraph'.format(str(window)))