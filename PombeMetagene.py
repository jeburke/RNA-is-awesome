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
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
import SPTools as SP
import PombePeaks as PP
import GeneTools as GT
from scipy.stats import ks_2samp
import itertools
import sklearn
from matplotlib import colors
from copy import deepcopy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import random
from scipy.cluster.hierarchy import leaves_list, linkage

### Transcript dictionary where key is transcript name and values are [start, stop, strand, chromosome, list of cds starts, list  of cds ends]. Can provide a different gff3 file if desired.
def build_transcript_dict(gff3="/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3", expand=False, convert_chroms=False):
    transcript_dict = SP.build_transcript_dict(gff3, organism='pombe')
    
    lat_rom = {'chr1':'I','chr2':'II','chr3':'III','MT':'MT'}
    
    if convert_chroms is True:
        transcript_dict = {k:[start, end, strand, lat_rom[chrom], cds_start, cds_end] for 
                           k, [start, end, strand, chrom, cds_start, cds_end] in transcript_dict.items()}
    
    
    chrom_lengths = {'I':5818680, 'II':4744158, 'III':2598968,'chr1':5818680, 'chr2':4744158, 'chr3':2598968}
    
    if expand is True:
        expanded_dict = {}
        for tx, info in transcript_dict.iteritems():
            new_start = info[0]-300
            if new_start < 0:
                new_start = 0
            new_end = info[1]+300
            if info[3] in chrom_lengths:
                if new_end > chrom_lengths[info[3]]:
                    new_end = chrom_lengths[info[3]]
            #else: print info[3]
            if len(info[4]) == 0:
                info[4] = [info[0]]
            if len(info[5]) == 0:
                info[5] = [info[1]]
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
    
## ****New code for travelling ratios starts here **** ##

def region_density(s, info, window_size=500):
    start, end, strand, CDS_start, CDS_end, exons, chrom = info
    start = start+300
    end = end-300
    ORF = s[s.index.isin(range(start,end))]
    if sum(ORF) > 50:
        if end-start >= 1000:
            if strand == '+':
                range1 = range(start,start+window_size)
                range2 = range(end-window_size,end)
                range3 = range(end,end+100)
            
            elif strand == '-':
                range1 = range(end-window_size,end)
                range2 = range(start,start+window_size)
                range3 = range(start-100,start)
            
            sums = {}
            ranges = {'five':(range1,float(window_size)), 'three':(range2,float(window_size)), 'ds':(range3,100.)}
            for name, (rg, window) in ranges.iteritems():
                rg_sum = sum(s[s.index.isin(rg)])
                if rg_sum == 0: rg_sum += 1
                sums[name] = rg_sum/float(window)
                
            denom = (float(sum(ORF))/(end-start))
            five = sums['five']/denom
            three = sums['three']/denom
            ds = sums['ds']/denom
        else: 
            return None
    else: 
        return None
        
    return five, three, ds

def heatmap_bins(s, info, window_size=500):
    start, end, strand, CDS_start, CDS_end, exons, chrom = info
    ORF = s[s.index.isin(range(start,end))]
    ORF_sum = float(sum(ORF))/(end-start)
    if sum(ORF) > 0:
        if end-start >= 1000:
            five_sums = []
            three_sums = []
            if strand == '+':
                range1 = range(start,start+window_size)
                range2 = range(end-window_size,end)
                
                n = min(range1)
                while n < max(range1):
                    bin_sum = sum(s[s.index.isin(range(n,n+100))])
                    if bin_sum == 0: bin_sum += 1
                    bin_sum = bin_sum/100.
                    five_sums.append(bin_sum/ORF_sum)
                    n += 100
                n = min(range2)
                while n < max(range2):
                    bin_sum = sum(s[s.index.isin(range(n,n+100))])
                    if bin_sum == 0: bin_sum += 1
                    bin_sum = bin_sum/100.
                    three_sums.append(bin_sum/ORF_sum)
                    n += 100
                
            elif strand == '-':
                range1 = range(end-window_size,end)
                range2 = range(start,start+window_size)
                
                n = max(range1)
                while n > min(range1):
                    bin_sum = sum(s[s.index.isin(range(n-100,n))])
                    if bin_sum == 0: bin_sum += 1
                    bin_sum = bin_sum/100.
                    five_sums.append(bin_sum/ORF_sum)
                    n = n-100
                n = max(range2)
                while n > min(range2):
                    bin_sum = sum(s[s.index.isin(range(n-100,n))])
                    if bin_sum == 0: bin_sum += 1
                    bin_sum = bin_sum/100.
                    three_sums.append(bin_sum/ORF_sum)
                    n = n-100
        else: 
            return None
    else: 
        return None
        
    return five_sums, three_sums


def add_cdf_to_plot(ax, value_lists, label_list, color_list, ks_list, log2=False):
    all_cdfs = []
    all_lists = []
    n = 0 
    
    for n in range(len(value_lists)):
        if log2 is True:
            new_list = [np.log2(x) for x in value_lists[n]]
        else: new_list = value_lists[n]
        new_list = [x for x in new_list if (str(x) != 'inf' and str(x) != '-inf' and str(x) != 'nan') ]
        all_lists = all_lists+new_list
        cumulative, base = SP.cdf_values(new_list)
        ax.plot(base[1:], cumulative, c=color_list[n], linewidth=3.0, label=label_list[n])
        all_cdfs.append(cumulative)
        
    xmin = np.percentile(all_lists, 1)
    xmax = np.percentile(all_lists, 99)
    ax.set_xlim([xmin,xmax])
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    
    if ks_list is not None:
        text = "p-values:    \n"+ks_list[0]+'    \n'+ks_list[1]+'    '
        if len(ks_list) == 4:
            text = text+'    \n'+ks_list[2]+'    \n'+ks_list[3]+'    '
        ax.annotate(text,  xy=(xmax,0.0), horizontalalignment='right', fontsize=12)
    
    return ax
    

def plot_traveling_ratio(bam_files, tx_dict=None, label_list=None, color_list=None, stranded=True, xlabel='traveling ratio', ylabel='Fraction of transcripts', ks_test=False, testing=False, both=True):
    
    if tx_dict is None:
        tx_dict = build_transcript_dict(expand=True)

        if testing is True:
            tx_iter = itertools.islice(tx_dict.items(), 0, 20)
            mini_tx_dict = {}
            for key, value in tx_iter:
                mini_tx_dict[key] = value
            tx_dict = remove_overlapping_transcripts(mini_tx_dict)
        
        else:
            tx_dict = remove_overlapping_transcripts(tx_dict)

    # First get read series for all transcripts in all samples 
    read_dict = {}
    TR_dict = {} #bam file, transcript then (5p, 3p, ds)
    hm_dict = {} #for heatmpa, bam file, transcript them (5p_bins, 3p_bins)
    
    print "Calculating traveling ratios..."
    for bam in bam_files:
        print bam
        read_dict[bam] = {}
        TR_dict[bam] = {}
        hm_dict[bam] = {}
        bam_reader = pysam.Samfile(bam)
        lat_rom = {'chr1':'I','chr2':'II','chr3':'III'}
        for tx, info in tx_dict.iteritems():
            if info[3] in lat_rom: chrom = lat_rom[info[3]]
            else: chrom = info[3]

            iterator = bam_reader.fetch(chrom, info[0], info[1])
            s = PP.generate_read_series(iterator, chrom, info[0], info[1], info[2])
            
            if stranded is False:
                switch = {'+':'-','-':'+'}
                iterator = bam_reader.fetch(chrom, info[0], info[1])
                s2 = PP.generate_read_series(iterator, chrom, info[0], info[1], switch[info[2]])
                s = s.add(s2, fill_value=0)
            
            read_dict[bam][tx] = s

            # Calculate traveling ratios
            info = PP.tx_info(tx, tx_dict)
            TR_values = region_density(s, info, window_size=500)

            if TR_values is not None:
                TR_dict[bam][tx] = TR_values
                
            # Get bins for heatmap
            hm_dict[bam][tx] = heatmap_bins(s, info, window_size=500)
                
    # Convert to lists for plotting
    names = []
    five_lists = []
    three_lists = []
    ds_lists = []
    for n, bam in enumerate(bam_files):
        names.append(bam.split('/')[-1].split('_sorted.bam')[0])
        five_lists.append([])
        three_lists.append([])
        ds_lists.append([])
        
        for tx, TRs in TR_dict[bam].iteritems():
            five_lists[n].append(TRs[0])
            three_lists[n].append(TRs[1])
            ds_lists[n].append(TRs[2])
                
    if color_list is None:
        color_list = ['0.3','0.6','orangered','orange','teal','turquoise']

    ks_lists = [None,None,None]
    if ks_test is True:
        for m, lists in enumerate([five_lists, three_lists, ds_lists]):
            ks_lists[m] = []
            ks_lists[m].append("%0.1e" % ks_2samp(lists[0], lists[2])[1])
            ks_lists[m].append("%0.1e" % ks_2samp(lists[1], lists[3])[1])
            if len(lists) == 6:
                ks_lists[m].append("%0.1e" % ks_2samp(lists[0], lists[4])[1])
                ks_lists[m].append("%0.1e" % ks_2samp(lists[1], lists[5])[1])
                
    wind_lists = []
    if both is False:
        wind_lists.append([five_lists[0],five_lists[2]])
        wind_lists.append([three_lists[0],three_lists[2]])
        wind_lists.append([ds_lists[0],ds_lists[2]])
        names = [names[0],names[2]]
    else:
        wind_lists = [five_lists, three_lists, ds_lists]

    # Set up figure and plot
    fig, (ax1, ax2) = plt.subplots(2, 3, figsize=(14, 10), dpi=600)
    
    for n, wind in enumerate(["5'", "3'","Downstream"]):
        ax1[n] = add_cdf_to_plot(ax1[n], wind_lists[n], names, color_list, ks_lists[n])
        ax1[n].set_ylabel(ylabel, fontsize=14)
        ax1[n].set_xlabel(wind+' '+xlabel+'\n500 nt window', fontsize=14)
        
        ax2[n] = add_cdf_to_plot(ax2[n], wind_lists[n], names, color_list, ks_lists[n], log2=True)
        ax2[n].set_ylabel(ylabel, fontsize=14)
        ax2[n].set_xlabel(wind+' log2 '+xlabel+'\n500 nt window', fontsize=14)
    
    # Draw legend after 3rd plot
    ax1[0].legend(fontsize=12)
    fig.tight_layout()
    plt.show()
    
    return fig, TR_dict, hm_dict

def compare_TR(wt1, wt2, mut1, mut2):
    up_downs_both = {"Decreased 5'":set(),"Increased 5'":set(),"Decreased 3'":set(),"Increased 3'":set()}
    for tx in wt1:
        try:
            ratio51 = np.log2(mut1[tx][0]/wt1[tx][0])
            ratio52 = np.log2(mut2[tx][0]/wt2[tx][0])
            if ratio51 >= 0.58 and ratio52 >= 0.58:
                up_downs_both["Increased 5'"].add(tx)
            elif ratio51 <= -0.58 and ratio52 <= -0.58:
                up_downs_both["Decreased 5'"].add(tx)
        
            ratio31 = np.log2(mut1[tx][1]/wt1[tx][1])
            ratio32 = np.log2(mut2[tx][1]/wt2[tx][1])
            if ratio31 >= 0.58 and ratio32 >= 0.58:
                up_downs_both["Increased 3'"].add(tx)
            elif ratio31 <= -0.58 and ratio32 <= -0.58:
                up_downs_both["Decreased 3'"].add(tx)

        except ZeroDivisionError:
            pass
        except KeyError:
            pass
    
    all_other = set(wt1.keys())
    for category, tx_set in up_downs_both.iteritems():
        all_other = all_other.difference(tx_set)
    print len(all_other)
    up_downs_both["Unchanged"] = all_other
    return up_downs_both
    

def TR_venns(bam_list, TR_dict, N=4675):
    '''
    Parameters
    ----------
    bam_list : list, should be [wt1, wt2, mut1, mut2]
    TR_dict : dictionary, traveling ratios from plot_traveling_ratio function
    N : number of genes in population (default excludes overalaping genes and ncRNAs)
    
    Returns
    -------
    up_downs_both : '''
    
    up_downs_both = compare_TR(TR_dict[bam_list[0]], TR_dict[bam_list[1]], TR_dict[bam_list[2]], TR_dict[bam_list[3]])
    
    combos = [("Increased 5'","Increased 3'"),("Decreased 5'","Decreased 3'"),
              ("Increased 5'","Decreased 3'"),("Decreased 5'","Increased 3'"),
             ("Increased 5'","Unchanged")]
    
    #print "*******Replicate 1 plots*******\n"
    classes = {}
    for a,b in combos:
        K = len(up_downs_both[a])
        J = len(up_downs_both[b])
        overlap = set(up_downs_both[a]).intersection(up_downs_both[b])
        k = len(overlap)
        n = K + J - k
        p = GT.hypergeometric(N,n,K,J,k)
        GT.venn_2sample(n, K, k, J, a, b, ['crimson','deepskyblue','darkorchid'], p)
        
        if a not in classes:
            classes[a] = up_downs_both[a]
        if b not in classes:
            classes[b] = up_downs_both[b]
        
        if p < 0.05:
            classes[a+' '+b] = overlap
            classes[a] = classes[a].difference(overlap)
            classes[b] = classes[b].difference(overlap)
    
    return classes

def TR_heatmap2(TR_dict, bam_files, clusters=None, name='HeatMap'):
    if clusters is not None:
        index = set()
        for cluster in clusters:
            index.update(cluster)
    else:
        index = TR_dict[bam_files[0]].keys()
    
    columns = pd.MultiIndex.from_product([['TSS','PolyA'],['WT1','WT2','Mut1','Mut2']])
    TR_df = pd.DataFrame(index=index, columns=columns)
    
    labels = ['WT1','WT2','Mut1','Mut2']
    for n,bam in enumerate(bam_files):
        for tx, TRs in TR_dict[bam].iteritems():
            if TRs is not None:
                five, three, ds = TRs
                TR_df.loc[tx, ('TSS',labels[n])] = five
                TR_df.loc[tx, ('PolyA',labels[n])] = three
    TR_df = TR_df.replace([np.inf,np.inf*-1],np.NaN)
    TR_df = TR_df.dropna(how='any')
    
    all_data = []
    for column in TR_df:
        all_data = all_data + TR_df[column].tolist()
    
    # Sort by clusters 
    cluster_sizes = []
    if clusters is not None:
        new_TR = None
        cluster_names = sorted(clusters.keys(), reverse=True)
        for name in cluster_names:
            cluster = clusters[name]
            cluster_df = TR_df[TR_df.index.isin(cluster)]
            if "5'" in name:
                cluster_df.loc[:,('TSS','ratio')] = ((cluster_df[('TSS','Mut1')]/cluster_df[('TSS','WT1')])+
                                                   (cluster_df[('TSS','Mut2')]/cluster_df[('TSS','WT2')])).divide(2.)
                cluster_df = cluster_df.sort_values(('TSS','ratio'), ascending=False)
                cluster_df = cluster_df.drop(('TSS','ratio'), axis=1)
            elif "3'" in name:
                cluster_df.loc[:,('PolyA','ratio')] = ((cluster_df[('PolyA','Mut1')]/cluster_df[('PolyA','WT1')])+
                                                   (cluster_df[('PolyA','Mut2')]/cluster_df[('PolyA','WT2')])).divide(2.)
                cluster_df = cluster_df.sort_values(('PolyA','ratio'), ascending=False)
                cluster_df = cluster_df.drop(('PolyA','ratio'), axis=1)
            elif name == 'Unchanged':
                if len(cluster_df) > 100:
                    random_index = random.sample(cluster_df.index, 100)
                    cluster_df = cluster_df[cluster_df.index.isin(random_index)]
                    cluster_df = cluster_df.sort_values(('TSS','WT1'))
                    
            cluster_sizes.append(len(cluster_df))
            
            if new_TR is None:
                new_TR = deepcopy(cluster_df)
            else:
                new_TR = pd.concat([new_TR, cluster_df])
    else:
        TR_df[('TSS','Ratio 1')] = TR_df[('TSS','Mut1')]/TR_df[('TSS','WT1')]
        TR_df[('PolyA','Ratio 1')] = TR_df[('PolyA','Mut1')]/TR_df[('PolyA','WT1')]
        TR_df[('TSS','Ratio 2')] = TR_df[('TSS','Mut2')]/TR_df[('TSS','WT2')]
        TR_df[('PolyA','Ratio 2')] = TR_df[('PolyA','Mut2')]/TR_df[('PolyA','WT2')]
        TR_df = TR_df.replace([np.inf,np.inf*-1],np.NaN).dropna(how='any')
        
        #mat = TR_df.as_matrix(columns=[('TSS','Ratio 2'),('PolyA','Ratio 2')])
        mat = TR_df.as_matrix(columns=[('TSS','WT2'),('TSS','Mut2'),('PolyA','WT2'),('PolyA','Mut2')])
        km = sklearn.cluster.KMeans(n_clusters=3, random_state=0).fit(mat)
        labels = km.labels_
        
        #Z = linkage(mat, 'ward')
        #print leaves_list(Z)
        
        TR_df[('All','cluster')] = labels
        cluster_names = list(set(labels))
        TR_df[('TSS','WT avg')] = (TR_df[('TSS','WT1')] + TR_df[('TSS','WT2')]).divide(2.)
        TR_df = TR_df.sort_values([('All','cluster'),('TSS','WT avg')])
        
        cluster_sizes = []
        for index in range(max(TR_df[('All','cluster')])):
            cluster_sizes.append(len(TR_df[TR_df[('All','cluster')] == index]))
            print index
            print len(TR_df[TR_df[('All','cluster')] == index])
            for column in TR_df.columns:
                print column
                print np.median(TR_df[TR_df[('All','cluster')] == index][column].apply(np.log2))
            
        cluster_keys = TR_df[('All','cluster')]
        TR_df = TR_df.drop([('All','cluster'),('TSS','Ratio 1'),('TSS','Ratio 2'),('PolyA','Ratio 1'),('PolyA','Ratio 2'),('TSS','WT avg')], axis=1)
        new_TR = deepcopy(TR_df)
    
    # Figure out scales
    data_min = np.percentile(all_data, 2)
    data_max = np.percentile(all_data, 98)
    both_max = max([data_min*-1, data_max])
    print both_max
    
    for column in new_TR.columns:
        new_TR[column] = pd.to_numeric(new_TR[column].apply(np.log2))
    
    # Make heatmap
    fig = plt.figure(figsize=(4,12))
    ax = fig.add_subplot(111)
    pcm = ax.pcolor(range(len(new_TR.columns)+1), range(len(new_TR.index)), new_TR, cmap='RdBu_r', 
                    norm=colors.Normalize(vmin=both_max*-1, vmax=both_max))

    tick_spacing = []
    base = 0
    for size in cluster_sizes:
        tick_spacing.append(base+size)
        base += size

    ax.yaxis.set_ticks(tick_spacing)
    ax.yaxis.set_ticklabels(cluster_names, fontsize=14)
    ax.set_xticklabels([])
    
    # Label replicates
    for n, label in {1:'WT',3:'seb1-1',5:'WT',7:'seb1-1'}.iteritems():
        ax.text(n, len(new_TR)+5, label, horizontalalignment='center', fontsize=14)
    ax.text(2, len(new_TR)+20, 'TSS', horizontalalignment='center', fontsize=16)
    ax.text(6, len(new_TR)+20, 'PolyA', horizontalalignment='center', fontsize=16)
    
    for size in tick_spacing:
        ax.plot([0,len(new_TR.columns)], [size,size], '-', color='0.5')
    ax.plot([4,4], [0,len(new_TR)], '-', color='0.5')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.3)
    cbar = plt.colorbar(pcm, cax=cax)
    cbar.ax.tick_params(labelsize=14)
    ax.text(8.2, len(new_TR)/2, 'log2 Traveling ratio', verticalalignment='center',fontsize=14,rotation='vertical')
    fig.tight_layout()
    fig.savefig(name+'.eps', format='eps')
    plt.show()
    
    if clusters is None:
        return fig, TR_df, cluster_keys
    else:
        return fig, TR_df

def unlog2(x):
    x = 2**x
    return x

