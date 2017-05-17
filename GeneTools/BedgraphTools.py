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

def generate_scaled_bedgraphs(directory, organism='crypto', start_only=False):
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
	aligned_reads = check_output(command.split(), shell=False)
	total_aligned.append(aligned_reads)
	print "Total aligned reads in "+bam
	print aligned_reads
    
    ratio_list = []
    n=0
    for n in range(len(bam_list)):
        print bam_list[n]
        print total_aligned[n]
        ratio = float(total_aligned[n])/float(total_aligned[0])
        ratio_list.append(1/ratio)
        
    for n in range(len(bam_list)):
        out = bam_list[n].split('/')[-1].split('.')[0]
        if start_only is False:
            command1 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand + -scale {2}'.format(bam_list[n], genome, str(ratio_list[n]))
        elif start_only is True:
            command1 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand + -5 -scale {2}'.format(bam_list[n], genome, str(ratio_list[n]))
        print command1
        bg1 = check_output(command1.split(), shell=False)
        with open('{0}_plus.bedgraph'.format(out),'w') as fout:
            fout.write(bg1)
        if start_only is False:
            command2 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand - -scale {2}'.format(bam_list[n], genome, str(ratio_list[n]))
        elif start_only is True:
            command2 = 'genomeCoverageBed -ibam {0} -g {1} -bg -strand - -5 -scale {2}'.format(bam_list[n], genome, str(ratio_list[n]))
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
