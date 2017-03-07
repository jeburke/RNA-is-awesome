import sys
import os
from subprocess import check_output
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import pysam
from scipy import stats
from statsmodels import robust
from scipy.interpolate import spline
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/')
import SPTools as SP

### Transcript dictionary where key is transcript name and values are [start, stop, strand, chromosome, list of cds starts, list  of cds ends]. Can provide a different gff3 file if desired.
def build_transcript_dict(gff3="/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3", expand=False):
    transcript_dict = SP.build_transcript_dict(gff3, organism='pombe')
    
    if expand is True:
        expanded_dict = {}
        for tx, info in transcript_dict.iteritems():
            new_start = info[0]-220
            new_end = info[1]+220
            expanded_dict[tx] = [new_start, new_end, info[2], info[3], info[4], info[5]]
        transcript_dict = expanded_dict
    
    return transcript_dict

def find_polyA_sites(transcript_dict, window=220):
    polyA_bg = SP.read_CNAGsort_bedgraph2('/home/jordan/GENOMES/POMBE/polyA_sites_CNAGsort.bedgraph', transcript_dict, organism='pombe')
    pA_dict = {}
    for tx, s in polyA_bg.iteritems():
        s = s[s > 0]
        if len(s) > 0:
            if transcript_dict[tx][2] == '+':
                #pA_site = max(s.index)
                s.sort_values(ascending=False, inplace=True)
                pA_site = s.index[0]
                pA_dict[tx] = [pA_site-window, pA_site+window, transcript_dict[tx][2], transcript_dict[tx][3]]
            elif transcript_dict[tx][2] == '-':
                #pA_site = min(s.index)
                s.sort_values(ascending=False, inplace=True)
                pA_site = s.index[0]
                pA_dict[tx] = [pA_site-window, pA_site+window, transcript_dict[tx][2], transcript_dict[tx][3]]
    return pA_dict

def build_tss_dict(gff3="/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3", window=220):
    transcript_dict = SP.build_transcript_dict(gff3, organism='pombe')
    
    tss_dict = {}
    for tx, info in transcript_dict.iteritems():
        if info[2] == '+':
            start = info[0]-window
            end = info[0]+window
            tss_dict[tx] = [start, end, info[2], info[3]]
        elif info[2] == '-':
            start = info[1]-window
            end = info[1]+window
            tss_dict[tx] = [start, end, info[2], info[3]]
    return tss_dict

def count_reads_in_window(bam, chrom, start, end, strand):
    bam_iter = bam.fetch(chrom, start, end)
    read_count = 0
    for read in bam_iter:
        if strand == "+":
            if read.is_reverse:
                read_count += 1
        elif strand == "-":
            if not read.is_reverse:
                read_count +=1
    return read_count

def count_aligned_reads(bam_file):
    bam = pysam.Samfile(bam_file)
    total = check_output(['samtools','view','-F 0x04','-c',bam_file]).strip()
    total = float(total)/1000000.
    return total

def rpkm_from_bam(transcript_dict, bam_file, fiveprime_only=False):
    bam = pysam.Samfile(bam_file)
    total = count_aligned_reads(bam_file)
    rpkm_tuples = []
    chrom_convert = {'chr1':'I','chr2':'II','chr3':'III'}
    
    for tx, info in transcript_dict.iteritems():
        chrom = info[3]
        if chrom in chrom_convert:
            chrom = chrom_convert[chrom]
        strand = info[2]
        start = info[0]
        if start < 0: start = 0
        end = info[1]
        if fiveprime_only is True:
            if strand == '+':
                end = start+200
            elif strand == '-':
                start = end-200
        length = abs(end-start)
        read_count = count_reads_in_window(bam, chrom, start, end, strand)
        if length >= 200:
            rpkm = float(read_count)/(length*total)
            rpkm_tuples.append((tx,rpkm))
    return rpkm_tuples

def find_abundant(rpkm_tuples, transcript_dict, cutoff=10):
    abundant = [x for x,y in rpkm_tuples if y > cutoff]
    abundant_dict = {k:v for k, v in transcript_dict.items() if k in abundant}
    return abundant_dict

#Need transcript list that is only non-overlapping genes with certain abundance
def remove_overlapping_transcripts(transcript_dict):
    no_overlap = transcript_dict
    for tx in transcript_dict.keys():
        if tx in no_overlap:
            for comp_tx in no_overlap.keys():
                if transcript_dict[tx][3] == no_overlap[comp_tx][3]:
                    if transcript_dict[tx][0] > no_overlap[comp_tx][0] and transcript_dict[tx][0] < no_overlap[comp_tx][1]:
                        deleted = no_overlap.pop(comp_tx, None)
                    elif transcript_dict[tx][1] >no_overlap[comp_tx][0] and transcript_dict[tx][1] < no_overlap[comp_tx][1]:
                        deleted = no_overlap.pop(comp_tx, None)
    print len(no_overlap)
    return no_overlap

def refine_transcript_dict(transcript_dict, rpkm_tuples, cutoff=10):
    abundant = find_abundant(rpkm_tuples, transcript_dict, cutoff=cutoff)
    no_overlap = remove_overlapping_transcripts(abundant)
    return no_overlap

def sort_bedgraphs(directory, transcript_dict):
    bedgraph_list = []
    for file in os.listdir(directory):
        if file.lower().endswith(".bedgraph"):
            print file
            bedgraph_list.append(directory+file)
            
    for bedgraph in bedgraph_list:
        SP.build_bedgraph_dict(transcript_dict, bedgraph)
        
def read_sorted_bedgraphs(directory, transcript_dict):
    stranded_bedgraphs = {}
    for file in os.listdir(directory):
        if file.endswith("_CNAGsort.bedgraph"):
            if "plus" in file:
                if file.split('_plus')[0] not in stranded_bedgraphs:
                    stranded_bedgraphs[file.split('_plus')[0]] = [None, None]
                stranded_bedgraphs[file.split('_plus')[0]][0] = SP.read_CNAGsort_bedgraph2(file, transcript_dict, organism='pombe')
            elif 'minus' in file:
                if file.split('_minus')[0] not in stranded_bedgraphs:
                    stranded_bedgraphs[file.split('_minus')[0]] = [None, None]
                stranded_bedgraphs[file.split('_minus')[0]][1] = SP.read_CNAGsort_bedgraph2(file, transcript_dict, organism='pombe')
    return stranded_bedgraphs

def build_metagene(tx_dict, bedgraph_dict_tuple, bam_file, window=440):
    metagene = np.zeros([1,window+1])[0]
    #aligned_reads_million = count_aligned_reads(bam_file)
    aligned_reads_million = 1
    
    for tx, data in bedgraph_dict_tuple[0].iteritems():
        if tx in tx_dict and tx_dict[tx][2] == '+':
            total = sum(data)
            start = tx_dict[tx][0]
            end = tx_dict[tx][1]
            middle = (end-start)/2+start
            region1 = data.ix[start:middle]
            region2 = data.ix[middle:end]
            n = 0
            for index, value in region1.iteritems():
                metagene[n] += float(value)/aligned_reads_million/total
                n += 1
            for index, value in region2.iteritems():
                try:
                    metagene[n] += float(value)/aligned_reads_million/total
                    n+=1
                except IndexError:
                    pass
   
    for tx, data in bedgraph_dict_tuple[1].iteritems():
        if tx in tx_dict and tx_dict[tx][2] == '-':
            data = data.sort_index(ascending=False)
            total = sum(data)
            start = tx_dict[tx][1]
            end = tx_dict[tx][0]
            middle = (start-end)/2+end
            region1 = data[start:middle]
            region2 = data[middle:end]
            n = 0
            for index, value in region1.iteritems():
                metagene[n] += float(value)/aligned_reads_million/total
                n += 1
            for index, value in region2.iteritems():
                try:
                    metagene[n] += float(value)/aligned_reads_million/total
                    n+=1
                except IndexError:
                    pass
    x = range(len(metagene))
    
    #x_smooth = np.linspace(min(x), max(x), 100)
    #y_smooth = spline(x, metagene, x_smooth)
    return (x, metagene)

def mean_confidence_interval(data, confidence=0.99):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), stats.sem(a)
    h = se * stats.t._ppf((1+confidence)/2., n-1)
    return m-h, m+h

def metagene_confidence_intervals(tx_dict, wt_bedgraph_dict_tuple, mut_bedgraph_dict_tuple, bam_file, window=600):
    wt_CI = np.zeros([2,window+220])
    
    n=0
    for n in range(window+219):
        position_list = []
        for tx, data in wt_bedgraph_dict_tuple[0].iteritems():
            if tx in tx_dict and tx_dict[tx][2] == '+':
                if len(data) == 820:
                    total = sum(data)
                    position = data.tolist()[n]/float(total)
                    position_list.append(position)
        
        for tx, data in wt_bedgraph_dict_tuple[1].iteritems():
            if tx in tx_dict and tx_dict[tx][2] == '-':
                if len(data) == 820:
                    total = sum(data)
                    data = data.sort_index(ascending=False)
                    position = data.tolist()[n]/float(total)
                    position_list.append(position)
        
        pos_CI = mean_confidence_interval(position_list)
        wt_CI[0,n] = pos_CI[0]
        wt_CI[1,n] = pos_CI[1]
    
    below_CI = []
    above_CI = []
    m=0
    
    for m in range(window+219):
        tx_count = 0
        position_sum = 0
        for tx, data in mut_bedgraph_dict_tuple[0].iteritems():
            if tx in tx_dict and tx_dict[tx][2] == '+':
                if len(data) == 820:
                    total = sum(data)
                    tx_count += 1
                    position = data.tolist()[m]/float(total)
                    position_sum += position
        
        for tx, data in mut_bedgraph_dict_tuple[1].iteritems():
            if tx in tx_dict and tx_dict[tx][2] == '-':
                if len(data) == 820:
                    total = sum(data)
                    tx_count += 1
                    data = data.sort_index(ascending=False)
                    position = data.tolist()[m]/float(total)
                    position_sum += position
        position_avg = float(position_sum)/tx_count
        #print wt_CI[0,m]
        #print wt_CI[1,m]
        #print position_avg
        #print '\n'
        
        if position_avg < wt_CI[0,m]:
            below_CI.append(-window+m)
        elif position_avg > wt_CI[1,m]:
            above_CI.append(-window+m)
    return below_CI, above_CI                 
                    
def metagenes_from_bedgraph(stranded_bedgraphs, transcript_dict, original_directory='./', window=440):
    metagene_dict = {}
    for sample in stranded_bedgraphs:
        bam_file = original_directory+sample+'_primed_sorted.bam'
        metagene_dict[sample] = build_metagene(transcript_dict, stranded_bedgraphs[sample], bam_file, window=window)
    return metagene_dict

def plot_metagenes(sample_names, metagene_dict, CI_lists=None, color_list=None, plot_name='Yay', ymax=None, xlabel=None, title=None):
    if color_list is None:
        color_list = ['0.5','orange','navy','skyblue']
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ymax_list = []
    xmin_list = []
    xmax_list = []
    n=0
    for n in range(len(sample_names)):
        sample = sample_names[n]
        half = (len(metagene_dict[sample][1])+1)/2
        x = range(3-half,0+half)
        print len(x)
        y = metagene_dict[sample][1][2:]
        print len(y)
        x_smooth = np.linspace(min(x), max(x), 100)
        y_smooth = spline(x, y, x_smooth)
        ymax_list.append(max(y))
        xmin_list.append(min(x))
        xmax_list.append(max(x))
        ax.plot(x_smooth, y_smooth, color_list[n], label=sample)
    
    if ymax is None:
        ymax = max(ymax_list)+0.1*max(ymax_list)
    xmin = max(xmin_list)+3
    xmax = min(xmax_list)-3
    if CI_lists is not None:
        y_low = []
        for x in CI_lists[0]:
            y_low.append(ymax)
        y_high = []
        for x in CI_lists[1]:
            y_high.append(ymax)
        plt.bar(CI_lists[0], y_low, color='0.9', alpha=0.1, zorder=10)
        plt.bar(CI_lists[1], y_high, color='red', alpha=0.1, zorder=10)
        

    ax.plot([0, 0], [0, ymax], color='0.7', linestyle='--', linewidth=2)
    if xlabel is None: ax.set_xlabel('Distance from gene feature')
    else: ax.set_xlabel(xlabel)
    ax.set_ylabel('Total normalized reads')
    if title is not None:
        ax.set_title(title)
    plt.xlim([xmin,xmax])
    plt.ylim([0,ymax])
    plt.legend()
    plt.show()
    fig.savefig(plot_name+'.pdf', format='pdf')
        
def plot_metagene_tss(sample_names, metagene_dict, color_list=None, plot_name='Yay', window=220):
    if color_list is None:
        color_list = ['0.5','orange','navy','skyblue']
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ymax_list = []
    xmin_list = []
    xmax_list = []
    n=0
    for n in range(len(sample_names)):
        sample = sample_names[n]
        x = [x+(window-220) for x in metagene_dict[sample][0]][1:]
        y = metagene_dict[sample][1][1:]
        ymax_list.append(max(y))
        xmin_list.append(min(x))
        xmax_list.append(max(x))
        plt.plot(x, y, color_list[n], label=sample)
        
    ymax = max(ymax_list)+0.5
    xmin = max(xmin_list)+3
    xmax = min(xmax_list)-3
    plt.plot([0, 0], [0, ymax], color='0.7', linestyle='--', linewidth=2)
    plt.xlim([xmin,xmax])
    plt.ylim([0,ymax])
    plt.legend()
    plt.show()
    fig.savefig(plot_name+'.pdf', format='pdf')        