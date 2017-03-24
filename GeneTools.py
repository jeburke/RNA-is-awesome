import sys
import os
from subprocess import check_output
from subprocess import call
from subprocess import Popen
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneUtility
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/SP_ANALYSIS/')
import SPTools as SP
from collections import OrderedDict
import csv
from itertools import islice


#Transcript dictionary: keys are transcript, values are [start, end, strand, chromosome, CDS start, CDS end]

def write_transcript_fasta(transcript_dict, fasta_dict, prefix='transcripts', sense=True, spliced=False):
    seq_dict = {}
    for transcript, values in transcript_dict.iteritems():
        start = values[0]
        end = values[1]
        strand = values[2]
        chrom = values[3]
        CDS_start_list = values[4]
        CDS_end_list = values[5]
        
        if spliced is False:
            seq = fasta_dict[chrom][start-1:end]
            if strand == '-':
                seq = SP.reverse_complement(seq)

        elif spliced is True:
            seq = ''
            for n in range(len(CDS_start_list)):
                if strand == '+':
                    seq = seq+fasta_dict[chrom][CDS_start_list[n]-1:CDS_end_list[n]]
                elif strand == '-':
                    new_seq = fasta_dict[chrom][CDS_start_list[n]-1:CDS_end_list[n]]
                    new_seq = SP.reverse_complement(new_seq)
                    seq = seq+new_seq
        
        if sense is False:
            seq = SP.reverse_complement(seq)
            
        seq_dict[transcript] = seq
        
    with open('{}.fa'.format(prefix), 'w') as fout:
        for transcript, seq in seq_dict.iteritems():
            fout.write('>'+transcript+'\n')
            fout.write(seq+'\n')
    
    return seq_dict
            
def write_intergenic_fasta(transcript_dict, fasta_dict, bps_us=0, bps_ds=0, all_intergenic=True, prefix='intergenic_transcripts'):
    seq_dict = {}
    if all_intergenic is False:
        for transcript, values in transcript_dict.iteritems():
            start = values[0]
            end = values[1]
            strand = values[2]
            chrom = values[3]
            
            if bps_us > 0:
                if strand == '+':
                    seq_us_sense = fasta_dict[chrom][start-bps_us:start]
                elif strand == '-':
                    seq_us_sense = fasta_dict[chrom][end:end+bps_us]
                    seq_us_sense = SP.reverse_complement(seq_us_sense)
                seq_us_antisense = SP.reverse_complement(seq_us_sense)
                seq_dict[transcript+'_us_sense'] = seq_us_sense
                seq_dict[transcript+'_us_antisense'] = seq_us_antisense
            
            if bps_ds > 0:
                if strand == '+':
                    seq_ds_sense = fasta_dict[chrom][end:bps_ds+end]
                elif strand == '-':
                    seq_ds_sense = fasta_dict[chrom][start-bps_ds:start]
                    seq_ds_sense = SP.reverse_complement(seq_ds_sense)
                seq_ds_antisense = SP.reverse_complement(seq_ds_sense)
                seq_dict[transcript+'_ds_sense'] = seq_ds_sense
                seq_dict[transcript+'_ds_antisense'] = seq_ds_antisense
    
    elif all_intergenic is True:
        chroms = fasta_dict.keys()
        for chrom in chroms:
            chrom_transcripts = dict((k, transcript_dict[k]) for k in transcript_dict if transcript_dict[k][3] == chrom)
            chr_txs_df = pd.DataFrame.from_dict(chrom_transcripts, orient='index')
            chr_txs_df.sort_values([0], inplace=True)
            sorted_transcripts = chr_txs_df.index.tolist()

            n = 0
            for n in range(len(sorted_transcripts)-1):
                transcript = sorted_transcripts[n]
                next_transcript = sorted_transcripts[n+1]
                transcript_end = chr_txs_df[1][transcript]
                next_start = chr_txs_df[0][next_transcript]
                if next_start > transcript_end:
                    seq_plus = fasta_dict[chrom][transcript_end:next_start]
                    seq_dict[transcript+'_'+next_transcript+'_plus'] = seq_plus
                    seq_dict[transcript+'_'+next_transcript+'_minus'] = SP.reverse_complement(seq_plus)
                else:
                    print 'Overlapping transcripts:'
                    print transcript
                    print next_transcript
                    
    with open('{}.fa'.format(prefix), 'w') as fout:
        for transcript, seq in seq_dict.iteritems():
            fout.write('>'+transcript+'\n')
            fout.write(seq+'\n')
    
    return seq_dict
        
def write_intron_fasta(transcript_dict, fasta_dict, prefix='introns', sense=True):
    seq_dict = {}
    for transcript, values in transcript_dict.iteritems():
        start = values[0]
        end = values[1]
        strand = values[2]
        chrom = values[3]
        CDS_start_list = values[4]
        CDS_end_list = values[5]

        for n in range(len(CDS_start_list)-1):
            if strand == '+':
                seq = fasta_dict[chrom][CDS_end_list[n]:CDS_start_list[n+1]-1]
            elif strand == '-':
                intron = len(CDS_start_list)-n-1
                seq = fasta_dict[chrom][CDS_end_list[intron]:CDS_start_list[intron-1]-1]
                seq = SP.reverse_complement(seq)
        
            if sense is False:
                seq = SP.reverse_complement(seq)
            
            seq_dict[transcript+'_'+str(n)] = seq
        
    with open('{}.fa'.format(prefix), 'w') as fout:
        for transcript, seq in seq_dict.iteritems():
            fout.write('>'+transcript+'\n')
            fout.write(seq+'\n')
    return seq_dict
        
    
def read_fasta(fasta_file):
    cds_dict = {}
    nts = ['A','T','C','G','N']
    with open(fasta_file, 'r') as fin:
        for line in fin:
            if line.startswith('>'):
                cds_dict[line.strip()] = ''
                cds = line.strip()
            elif line[0] in nts:
                cds_dict[cds] = cds_dict[cds]+line.strip()
    return cds_dict


def build_transcript_dict(gff3_file, organism=None):
    with open(gff3_file,"r") as gff3:
        transcript_dict = {}
        for line in gff3:
            columns = re.split(r'\t+', line)
            
            if organism == 'pombe' and len(columns) > 1:
                chr_rom = columns[0]
                rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
                chrom = rom_lat[chr_rom]
                transcript_types = ['transcript','pseudogene','rRNA','snoRNA','tRNA','snRNA']
                if columns[2] in transcript_types:
                    if columns[8].split(':')[0].split('=')[1] == 'gene': continue
                    transcript = columns[8].split(';')[0].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], chrom, [], []]
                elif columns[2] == 'CDS':
                    transcript = columns[8].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
                        
            if len(columns) == 9 and organism is None:
                if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                    transcript = columns[8]
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"): transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    #Transcript dictionary: keys are transcript, values are [start, end, strand, chromosome, CDS start, CDS end]
                    if transcript not in transcript_dict:
                        transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], columns[0], [], []]
                    else:
                        transcript_dict[transcript][0] = int(columns[3])
                        transcript_dict[transcript][1] = int(columns[4])
                        transcript_dict[transcript][2] = columns[6]
                        transcript_dict[transcript][3] = columns[0]
                elif columns[2] == "CDS":
                    transcript = columns[8].split("=")[1].split(".")[0].split(';')[0]
                    if 'mRNA' in transcript: transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    if transcript not in transcript_dict:
                        strand = columns[6]
                        chrom = columns[0]
                        transcript_dict[transcript] = [0,0,strand,chrom,[],[]]
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
    
    for tx in transcript_dict:
        if transcript_dict[tx][0] == 0:
            transcript_dict[tx][0] = transcript_dict[tx][4][0]
            transcript_dict[tx][1] = transcript_dict[tx][5][0]
    
    transcript_dict = collections.OrderedDict(sorted(transcript_dict.items()))
    return transcript_dict

def split_cds(cds_dict, end, number_nts):
    split_dict = {}
    for cds, seq in cds_dict.iteritems():
        if len(seq) <= number_nts:
            split_dict[cds] = seq
            continue
        if end == 'five_prime':
            split_dict[cds.split(' ')[0]+'_5'] = seq[:number_nts+1]
            split_dict[cds.split(' ')[0]+'_3'] = seq[number_nts+1:]
        elif end == 'three_prime':
            cds_len = len(seq)
            split_dict[cds.split(' ')[0]+'_5'] = seq[:(cds_len-number_nts)]
            split_dict[cds.split(' ')[0]+'_3'] = seq[(cds_len-number_nts):]
    return split_dict

def write_new_fasta(split_dict, fasta_file):
    new_name = fasta_file.split('.')[0]+'_split.fa'
    with open(new_name, 'w') as fout:
        for cds, seq in split_dict.iteritems():
            fout.write(cds+'\n')
            fout.write(seq+'\n')

def seq_simple(chrom, start, end, strand, fasta_dict):
    seq = fasta_dict[chrom][start:end+1]
    if strand == '-':
        seq = SP.reverse_complement(seq)
    return seq
            
def get_peak_sequence(input_file, fasta_file, gff3_file, window=1000):
    #File where 1st column is transcript (3P prefix sometimes), 2nd column is chr, 3rd column is peak center
    tx_dict = SP.build_transcript_dict(gff3_file)
    if type(fasta_file) == dict:
	fa_dict = fasta_file
    else:
        fa_dict = SP.make_fasta_dict(fasta_file)
    seq_list = []
    with open(input_file,'r') as csv_file:
        f = csv.reader(csv_file, dialect=csv.excel)
        for row in f:
            tx = row[0]+'T0'
            if tx.startswith('3P'): tx = tx.split('3P')[1]
            chrom = 'chr'+row[1]
            try:
                center = int(row[2])
                if tx in tx_dict:
                    strand = tx_dict[tx][2]
                    start = center-window/2
                    end = center+window/2
                    seq = seq_simple(chrom, start, end, strand, fa_dict)
                    seq_list.append((tx,seq))
                else:
                    print tx+" not in GFF3 file"
            except ValueError:
		pass
    with open('{0}_peak_sequences.fa'.format(input_file.split('/')[-1].split('.')[0]),'w') as fout:
        for tx, seq in seq_list:
	    fout.write('>'+tx+'\n')
	    fout.write(seq+'\n')
    return seq_list


# This function reads my sorted-by-gene bedgraph format (output of SPTools.build_bedgraph_dict). Make sure the same
# transcript_dictionary is used for both function calls. In this case I have modified the dictionary to include 300 bp
# downstrame of the stop codon for each transcript

def read_bg_dict(bedgraph_dict_output, transcript_dict):
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
                    if tx[-2] != 'T':
                        tx = tx+'T0'
                
                #Read the coordinate line
                elif n%3 == 2:
                    coords  = map(int, line.strip().split('\t'))
                
                #Read the values line
                elif n%3 == 0:
                    heights = map(float, line.strip().split('\t'))

                    all_coords = set(range(transcript_dict[tx][0],transcript_dict[tx][1]))
                    missing = all_coords.difference(coords)
                    coords = coords + list(missing)

                    #Fill in missing coordinates with zeros
                    zero_fill = [0]*len(missing)
                    heights = heights + zero_fill
                    
                    #Create a pandas series with all coordinates and sort so zeros are inserted appropriately
                    entry = pd.Series(heights, index=coords)
                    entry.sort_index(inplace=True)

                    bg_dict[tx] = entry
    return bg_dict

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

def generate_scaled_bedgraphs(directory, organism='crypto'):
    if 'crypto' in organism.lower():
        genome = '/home/jordan/GENOMES/crypto_for_bedgraph.genome'
    elif 'cerev' in organism.lower():
        genome = '/home/jordan/GENOMES/S288C/S288C_for_bedgraph.genome'
    elif 'pombe' in organism.lower():
        genome = '/home/jordan/GENOMES/POMBE/Sp_for_bg.genome'
    
    bam_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("sorted.bam"):
            bam_list.append(directory+file)
            
    total_aligned = []
    for bam in bam_list:
        command = 'samtools view -F 0x904 -c {0}'.format(bam)
        total_aligned.append(check_output(command.split(), shell=False))
    
    ratio_list = []
    n=0
    for n in range(len(bam_list)):
        print bam_list[n]
        print total_aligned[n]
        ratio = float(total_aligned[n])/float(total_aligned[0])
        ratio_list.append(1/ratio)
        
    for n in range(len(bam_list)):
        out = bam_list[n].split('/')[-1].split('.')[0]
        command1 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand + -scale {2}'.format(bam_list[n], genome, str(ratio_list[n]))
        print command1
        bg1 = check_output(command1.split(), shell=False)
        with open('{0}_plus.bedgraph'.format(out),'w') as fout:
            fout.write(bg1)
        command2 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand - -scale {2}'.format(bam_list[n], genome, str(ratio_list[n]))
        bg2 = check_output(command2.split(), shell=False)
        with open('{0}_minus.bedgraph'.format(out),'w') as fout:
            fout.write(bg2)

            
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