import sys
import os
from subprocess import check_output
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import pysam
from scipy import stats
from statsmodels import robust
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/')
import SPTools as SP
import random


### Transcript dictionary where key is transcript name and values are [start, stop, strand, chromosome, list of cds starts, list  of cds ends]. Can provide a different gff3 file if desired.
def build_transcript_dict(gff3="/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3"):
    transcript_dict = SP.build_transcript_dict(gff3, organism='pombe')
    return transcript_dict

### Same as transcript dict but with whole centromeres and subregions
def build_centromere_dict(gtf='/home/jordan/GENOMES/POMBE/Centromere_annotations.gtf'):
    cen_dict = {}
    with open(gtf,'r') as f:
        for line in f:
            data = line.split('\t')
            chrom = data[0]
            start = int(data[3])
            end = int(data[4])
            strand = data[6]
            ID = data[8].split('"')[1]
            if strand == '+':
                ID = ID+'_plus'
            elif strand == '-':
                ID = ID+'_minus'
            if ID not in cen_dict:
                cen_dict[ID] = [start, end, strand, chrom]
            else:
                new_ID = ID+':'+str(start)+'-'+str(end)
                cen_dict[new_ID] = [start, end, strand, chrom]
    print len(cen_dict)

    for ID in cen_dict.keys():
        if ID.endswith('plus'):
            if ID.split('plus')[0]+'minus' not in cen_dict.keys():
                cen_dict[ID.split('plus')[0]+'minus'] = [cen_dict[ID][0], cen_dict[ID][1], '-', cen_dict[ID][3]]
        elif ID.endswith('minus'):
            if ID.split('minus')[0]+'plus' not in cen_dict.keys():
                cen_dict[ID.split('minus')[0]+'plus'] = [cen_dict[ID][0], cen_dict[ID][1], '+', cen_dict[ID][3]]

    cen_dict['plus_cDNA_4'] = [3762489, 3763935, '+', 'chr1',[3762489], [3763935]]
    cen_dict['minus_cDNA_4'] = [3762489, 3763935, '-', 'chr1',  [3762489], [3763935]]
    cen_dict['CenI_L_plus'] = [3753687, 3764532, '+', 'chr1', [3753687], [3764532]]
    cen_dict['CenI_L_minus'] = [3753687, 3764532, '-', 'chr1', [3753687], [3764532]]
    cen_dict['CenI_R_plus'] = [3777597, 3789421, '+', 'chr1', [3777597], [3789421]]
    cen_dict['CenI_R_minus'] = [3777597, 3789421, '-', 'chr1', [3777597], [3789421]]

    cen_dict['CenII_L_plus'] = [1602619, 1618230, '+', 'chr2', [1602619], [1618230]]
    cen_dict['CenII_L_minus'] = [1602619, 1618230, '-', 'chr2', [1602619], [1618230]]
    cen_dict['CenII_R_plus'] = [1630357, 1643789, '+', 'chr2', [1630357], [1643789]]
    cen_dict['CenII_R_minus'] = [1630357, 1643789, '-', 'chr2', [1630357], [1643789]]

    cen_dict['CenIII_L_plus'] = [1071454, 1092169, '+', 'chr3', [1071454], [1092169]]
    cen_dict['CenIII_L_minus'] = [1071454, 1092169, '-', 'chr3', [1071454], [1092169]]
    cen_dict['CenIII_R_plus'] = [1106673, 1137549, '+', 'chr3', [1106673], [1137549]]
    cen_dict['CenIII_R_minus'] = [1106673, 1137549, '-', 'chr3', [1106673], [1137549]]

    print len(cen_dict)
    return cen_dict
    
#Function to grab all transcript information including exons from transcript dict
def tx_info(tx, tx_dict):
    strand = tx_dict[tx][2]
    start = tx_dict[tx][0]
    end = tx_dict[tx][1]
    exons = None
    if len(tx_dict[tx]) > 4 and len(tx_dict[tx][4]) > 0:
        CDS_start = min(tx_dict[tx][4])
        CDS_end = max(tx_dict[tx][5])
        exons = zip(tx_dict[tx][4],tx_dict[tx][5])
    else:
        CDS_start = start
        CDS_end = end
        
    return start, end, strand, CDS_start, CDS_end, exons
    
    
### Method for determining Z scores on somewhat skewed distributions
def robust_z(x, med, mad):
    z = abs(x-med)/mad
    return z

### Builds a pandas series from a bam file for a set of coordinates. Index is genome coordinate and value is number of reads
def generate_read_series(bam_iterator, chrom, start, end, strand, baseline=0):
    s = pd.Series(baseline, index=range(start-25, end+25))
    for read in bam_iterator:
        if read.is_reverse and strand == '+':
            pos = read.reference_start
            if pos not in s.index:
                s[pos] = 0
            s[pos] += 1
        elif not read.is_reverse and strand == '-':
            pos = read.reference_start
            if pos not in s.index:
                s[pos] = 0
            s[pos] = s[pos]+1
    s = s.dropna()
    s = s[s > 0]
    s = s.sort_index()
    return s

### Find peaks using a Z-score cutoff within a given region. Automatically sets minimum Z-score at 2.
def find_peaks(s, use_robust=True, min_z=2, use_log=True):   
    if use_log is True:
        new_s = s.apply(np.log)
    else:
        #new_s = s
        new_s = s[s <= s.quantile(0.95)]
    if use_robust is False:
        s_Z = pd.Series(stats.mstats.zscore(new_s), index=new_s.index)
    elif use_robust is True:
        s_MAD = robust.mad(new_s.tolist())
        s_med = np.median(new_s.tolist())
        s_Z = new_s.apply(robust_z, args=(s_med, s_MAD))
    
    marker_list = [[],[]]
    for pos, score in s_Z.iteritems():
        if score > min_z and s[pos] > 10:
        #if score > min_z:
            marker_list[0].append(pos)
            marker_list[1].append(s[pos]+0.01*s[pos])
    return marker_list

### Pick peaks in NET-seq data. Also consolidates clusters of peaks into a single center (must be within 10 bp). Returns a peak dictionary. Keys are transcripts (or regions) and values are a tuple. The first position is a read series for the region for plotting. The second position is a pair of lists that are the x and y values for each peak picked (+10% for easy plotting).
def NETseq_peaks(bam_file, transcript_dict, bam2=None, use_robust=True, min_z=2, use_log=True, scramble=False, collapse=False):
    peak_dict = {}
    bam_reader = pysam.Samfile(bam_file)
    chrom_convert = {'chr1':'I','chr2':'II','chr3':'III'}
    for tx, info in transcript_dict.iteritems():
        chrom = info[3]
        if chrom in chrom_convert:
            chrom = chrom_convert[chrom]
        start = info[0]
        end = info[1]
        strand = info[2]
        iterator = bam_reader.fetch(chrom, start, end)
        s = generate_read_series(iterator, chrom, start, end, strand)
        
        if scramble is True:
            s.index = random.sample(range(min(s.index),max(s.index)),len(s))
        
        marker_list = find_peaks(s, use_robust=use_robust, min_z=min_z, use_log=use_log)
        
        if collapse is True:
            if len(marker_list[0]) > 0:
                new_markers = [[],[]]
                prox = [[],[]]
                n=0
                while n < len(marker_list[0]):
                    if n == len(marker_list[0])-1:
                        if len(prox[0]) == 0:
                            new_markers[0].append(marker_list[0][n])
                            new_markers[1].append(marker_list[1][n])
                        else:
                            prox[0].append(marker_list[0][n])
                            prox[1].append(marker_list[1][n])
                            new_markers[0].append(sum(prox[0])/len(prox[0]))
                            new_markers[1].append(max(prox[1]))
                    else:
                        if marker_list[0][n+1]-marker_list[0][n] < 10:
                            prox[0].append(marker_list[0][n])
                            prox[1].append(marker_list[1][n])
                        else:
                            if len(prox[0]) == 0:
                                new_markers[0].append(marker_list[0][n])
                                new_markers[1].append(marker_list[1][n])
                            else:
                                prox[0].append(marker_list[0][n])
                                prox[1].append(marker_list[1][n])
                                new_markers[0].append(sum(prox[0])/len(prox[0]))
                                new_markers[1].append(max(prox[1]))
                            prox = [[],[]]
                    n += 1
                marker_list = new_markers
        peak_dict[tx] = (s,marker_list)
    return peak_dict


### Find all peaks in PARALYZER output. Returns a peak dictionary. Keys are transcripts or regions and values are read series
def PARCLIP_peaks(bam_file, transcript_dict):
    peak_dict = {}
    bam_reader = pysam.Samfile(bam_file)
    chrom_convert = {'chr1':'I','chr2':'II','chr3':'III'}
    for tx, info in transcript_dict.iteritems():
        chrom = info[3]
        if chrom in chrom_convert:
            chrom = chrom_convert[chrom]
        start = info[0]
        end = info[1]
        strand = info[2]
        
        #Need to switch strand here because PAR_Clip reads are on the first strand
        switch_strand = {'+':'-','-':'+'}
        strand = switch_strand[strand]
        iterator = bam_reader.fetch(chrom, start, end)
        s = generate_read_series(iterator, chrom, start, end, strand)
        s = s[s > 0]
        if len(s) > 0:
            peak_dict[tx] = s
    return peak_dict

# Compare PAR-CliP and NET-seq peaks. Look for PAR-ClIP peak within 60 upstream of any NET-seq peak
def find_PC_NS_overlap(PC_dict, NS_dict, transcript_dict):
    peaks_in_both = set(PC_dict.keys()).intersection(NS_dict.keys())
    print "Transcripts with peaks in PAR-ClIP and NET-seq: "+str(len(peaks_in_both))
    overlap = {}
    for tx in peaks_in_both:
        #Dictionary where values are (matched_peaks, all_NS_peaks)
        overlap[tx] = [[],[]]
        strand = transcript_dict[tx][2]
        
        for NS_peak in NS_dict[tx][1][0]:
            overlap[tx][1].append(NS_dict[tx][1][0])
            for PC_peak in PC_dict[tx].index:
                if strand == '+':
                    if NS_peak in range(PC_peak,PC_peak+60):
                        overlap[tx][0].append(NS_peak)
                        break
                elif strand == '-':
                    if NS_peak in range(PC_peak-60,PC_peak):
                        overlap[tx][0].append(NS_peak)
                        break
    return overlap


def compare_peak_overlap(wt_NS, mut_NS, PC_dict, transcript_dict):
    wt_peaks = find_PC_NS_overlap(PC_dict, wt_NS, transcript_dict)
    mut_peaks = {}
    for tx, peaks in wt_peaks.iteritems():
        mut_peaks[tx] = [[],[mut_NS[tx][1][0]]]
        for peak in peaks[0]:
            ext_peak = range(peak-50, peak+50)
            for pos in ext_peak:
                if pos in mut_NS[tx][1][0]:
                    mut_peaks[tx][0].append(peak)
    return mut_peaks
                
def plot_overlap(sample_names, region1_NS, region1_PC, transcript_dict, region2_NS=None, region2_PC=None, transcript_dict2=None, fig_name='PAR_NET_peak_overlap'):
    contingency = pd.DataFrame(index=sample_names, columns=['Matched','Unmatched'])
    x = np.arange(len(sample_names))
    y = []
    n=0
    for n in range(len(region1_NS)):
        print '---------------------------------------------------------------------------------------'
        print '\nCoding gene analysis for '+sample_names[n]+'\n'
        
        if n == 0:
            overlap = find_PC_NS_overlap(region1_PC, region1_NS[n], transcript_dict)
        else:
            overlap = compare_peak_overlap(region1_NS[0], region1_NS[n], region1_PC, transcript_dict)
    
        matched = [len(v[0]) for k,v in overlap.items() if v[0] > 0]
        #print "\nNumber of NET-seq peaks adjacent to a PAR-ClIP peak: "+str(sum(matched))
        all_NS = [len(v[1]) for k,v in overlap.items()]
        all_NS = sum(all_NS)
        unmatched = all_NS-sum(matched)
        #print "\nTotal number of NET-seq peaks: "+str(all_NS)
        #print str(float(sum(matched))/all_NS*100)+'%'
        y.append(float(sum(matched))/all_NS)
        contingency.loc[sample_names[n],'Matched'] = sum(matched)
        contingency.loc[sample_names[n],'Unmatched'] = unmatched
        
    print contingency
    fig, ax = plt.subplots()
    width=0.35
    rects1 = ax.bar(x, y, width, color='darkslateblue',label='Coding genes')
    autolabel(rects1, ax)
    chi2, p, dof, expected = stats.chi2_contingency(contingency)
    print "p-value for region1: "+str(p)
    
    
    if region2_NS is not None:
        contingency2 = pd.DataFrame(index=sample_names, columns=['Matched','Unmatched'])
        y2 = []
        n=0
        for n in range(len(region2_NS)):
            print '--------------------------------------------------------------------------------------'
            print '\nCentromere analysis for '+sample_names[n]+'\n'
            if n == 0:
                overlap = find_PC_NS_overlap(region2_PC, region2_NS[n], transcript_dict2)
            else:
                overlap = compare_peak_overlap(region2_NS[0], region2_NS[n], region2_PC, transcript_dict2)

            matched = [len(v[0]) for k,v in overlap.items() if v[0] > 0]
            #print "\nNumber of NET-seq peaks adjacent to a PAR-ClIP peak: "+str(sum(matched))
            all_NS = [len(v[1]) for k,v in overlap.items()]
            all_NS = sum(all_NS)
            unmatched = all_NS-sum(matched)
            #print "\nTotal number of NET-seq peaks: "+str(all_NS)
            #print str(float(sum(matched))/all_NS*100)+'%'
            y2.append(float(sum(matched))/all_NS)
            contingency2.loc[sample_names[n],'Matched'] = sum(matched)
            contingency2.loc[sample_names[n],'Unmatched'] = unmatched
            
        rects2 = ax.bar(width+x, y2, width, color='skyblue', label='Centromeres')
        ax.set_xticks(x + width/2)
        autolabel(rects2, ax)
        
        print contingency2
        chi2, p, dof, expected = stats.chi2_contingency(contingency2)
        print "p-value for region2: "+str(p)
    
    ax.set_ylim([0,1])
    ax.set_ylabel('Proportion NET-seq peaks with upstream PAR-ClIP peak')
    ax.set_xticklabels(sample_names)
    fig.savefig(fig_name+'.pdf', format='pdf')


def autolabel(rects, ax, sci_not=False):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        if sci_not is False:
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, '%.2f' % height, ha='center', va='bottom')
        else:
            ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, '%.2E' % height, ha='center', va='bottom')

def plot_NETvsPAR(region_dict, PC_dict, NS_dict, save_dir='.', filt=None):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    #Option to filter for just a subset of regions
    if filt is not None:
        region_dict = {k:v for k, v in region_dict.items() if k.startswith(filt)}
    
    for frag in region_dict:
        if sum(NS_dict[frag][0]) >= 50:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(311)
                
            norm_s = NS_dict[frag][0]
            ax.bar(norm_s.index, norm_s, width=30, color='0.3', label='NETseq')
            y_max = max(norm_s)+0.1*max(norm_s)
            
            if len(NS_dict[frag][1][0]) > 0:
                markers = NS_dict[frag][1][1]
                ax.plot(NS_dict[frag][1][0], markers,'*', markersize=10, color='red', label='Called NETseq peaks',zorder=2)
                y_max = max(markers)+0.1*max(markers)

            if frag in PC_dict:
                ax2 = plt.subplot(312, sharex=ax)
                ax2.bar(PC_dict[frag].index, PC_dict[frag], width=30, color='skyblue', edgecolor='skyblue', label='SEB1 PAR-ClIP', zorder=1)
                ax2.set_ylabel('PAR-ClIP reads')
            
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_ylim(0,y_max)
            ax.set_ylabel('NET-seq reads')
            plt.title(frag, y=2.2)
            plt.show()
            fig.savefig('{0}/{1}_NET_PAR.pdf'.format(save_dir,frag), format='pdf')
            plt.clf()

def plot_NETseq(region_dict, NS_dict, save_dir='.', NS_dict2=None, bam1=None, bam2=None, filt=None):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    if filt is not None:
        region_dict = {k:v for k, v in region_dict.items() if k.startswith(filt)}
    
    #If plotting 2 bam files, grab the total number of reads for normalization
    if NS_dict2 is None:
        mill_reads1 = 1
    else:
        mill_reads1 = check_output(['samtools','view','-F 0x04','-c',bam1]).strip()
        mill_reads1 = float(mill_reads1)/1000000.
        print mill_reads1
        mill_reads2 = check_output(['samtools','view','-F 0x04','-c',bam2]).strip()
        mill_reads2 = float(mill_reads2)/1000000.
        print mill_reads2
        
    for frag in region_dict:
        if sum(NS_dict[frag][0]) >= 50:
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(311)
                
            norm_s = NS_dict[frag][0]/mill_reads1
            ax.bar(norm_s.index, norm_s, width=30, color='0.3', label='NETseq')
            y_max = max(norm_s)+0.1*max(norm_s)
            
            #if len(NS_dict[frag][1][0]) > 0:
            #    markers = [x/mill_reads1 for x in NS_dict[frag][1][1]]
            #    ax.plot(NS_dict[frag][1][0], markers,'*', markersize=10, color='red', label='Called NETseq peaks',zorder=2)
            #    y_max = max(markers)+0.1*max(markers)
            
            if NS_dict2 is not None:
                ax_mut = plt.subplot(312, sharex=ax)
                norm_s2 = NS_dict2[frag][0]/mill_reads2
                ax_mut.bar(norm_s2.index, norm_s2, color='0.5', width=30, edgecolor='0.5', alpha=0.6,zorder=3)
                if y_max < max(norm_s2):
                    y_max = max(norm_s2)+0.1*max(norm_s2)
                ax_mut.set_ylim(0,y_max)
                ax_mut.set_ylabel('Mutant NET-seq reads')
            
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_ylim(0,y_max)
            ax.set_ylabel('NET-seq reads')
            plt.title(frag, y=2.2)
            plt.show()
            fig.savefig('{0}/{1}_NETseq.pdf'.format(save_dir, frag), format='pdf')
            plt.clf()
            
def net_seq_near_PAR(sample_names, bam_list, region1_NS, region1_PC, transcript_dict, region2_NS=None, region2_PC=None, transcript_dict2=None, fig_name='PAR_NET_peak_overlap'):
    mill_reads_list = []
    for m, bam in enumerate(bam_list):
        mill_reads = check_output(['samtools','view','-F 0x04','-c',bam]).strip()
        mill_reads = float(mill_reads)/1000000.
        print mill_reads
        mill_reads_list.append(mill_reads)

    y=[]
    n=0
    for n, NS in enumerate(region1_NS):
        y.append(0)
        for tx, s in region1_PC.iteritems():
            for index, value in s.iteritems():
                par_range = range(index-50, index+50)
                ns_signal = NS[tx][0]
                ns_in_range = sum(ns_signal[ns_signal.index.isin(par_range)])
                y[n] += float(ns_in_range)/mill_reads_list[n]
        print y[n]
    
    y =[v/y[0] for v in y]
    x = np.arange(len(sample_names))  
    fig, ax = plt.subplots()
    width=0.35
    rects1 = ax.bar(x, y, width, color='darkslateblue',label='Coding genes')
    #autolabel(rects1, ax, sci_not=True)
    
    if region2_NS is not None:
        y2 = []
        ns_near_par2 = {}
        n=0
        for n, NS in enumerate(region2_NS):
            y2.append(0)
            for tx, s in region2_PC.iteritems():
                for index, value in s.iteritems():
                    par_range = range(index-50, index+50)
                    ns_signal = NS[tx][0]
                    ns_in_range = sum(ns_signal[ns_signal.index.isin(par_range)])
                    y2[n] += float(ns_in_range)/mill_reads_list[n]
            print y2[n]
    
        y2 =[v/y2[0] for v in y2]
        rects2 = ax.bar(width+x, y2, width, color='skyblue', label='Centromeres')
        ax.set_xticks(x + width/2)
        #autolabel(rects2, ax, sci_not=True)
    
    #ax.set_ylim([0,1])
    ax.set_ylabel('NET-seq signal surrounding PAR-ClIP peak')
    ax.set_xticklabels(sample_names)
    fig.savefig(fig_name+'.pdf', format='pdf')


#Function for finding overlapping peaks in data sets. If from_reps is True, averages the data between the two sets    
def compare_NETseq_peak_counts(wt, mut, from_reps=False):
    total1 = 0
    for tx, info in wt.iteritems():
        total1 += sum(info[0])
        #print total1
    total2 = 0
    for tx, info in mut.iteritems():
        total2 += sum(info[0])
    
    total1 = total1/1000000.
    total2 = total2/1000000.
    
    in_both = {}
    both_count = 0
    only_WT = {}
    only_WT_count = 0
    only_mut = {}
    only_mut_count = 0
    for tx, info in wt.iteritems():
        peaks = info[1][0]
        heights = info[1][1]
        try:
            mut_peaks = mut[tx][1][0]
            #mut_heights = mut[tx][1][1]
        except KeyError:
            mut_peaks = []
            #mut_heights = []
            
        if from_reps is False:
            s = info[0]
            try:
                s2 = mut[tx][0]
                ratio = float(sum(s))/sum(s2)
                #print ratio
            except ZeroDivisionError:
                ratio = 0
            except KeyError:
                ratio = 0
        elif from_reps is True:
            s = (info[0]/total1+mut[tx][0]/total2)/2.
            s = s.fillna(0)
            s2 = s
            ratio = 1
        
        in_both[tx] = [s,([],[])]
        only_WT[tx] = [s,([],[])]
        only_mut[tx] = [s,([],[])]
        if len(peaks) > 0 and ratio <= 5 and ratio >= 0.2:
            for mut_peak in mut_peaks:
                #index = mut_peaks.index(mut_peak)
                if mut_peak not in peaks:
                    only_mut[tx][1][0].append(mut_peak) 
                    only_mut[tx][1][1].append(s2[mut_peak])
                    only_mut_count += 1
                else:
                    both_count += 1
                    in_both[tx][1][0].append(mut_peak) 
                    in_both[tx][1][1].append(s2[mut_peak])
                                
            for wt_peak in peaks:
                index = peaks.index(wt_peak)
                if wt_peak not in mut_peaks:    
                    only_WT[tx][1][0].append(wt_peak)
                    only_WT[tx][1][1].append(s[wt_peak])
                    only_WT_count += 1

    print "Peaks in both samples:"
    print both_count
    print "Peaks only in 1st sample:"
    print only_WT_count
    print "Peaks only in 2nd sample:"
    print only_mut_count
    return in_both, only_WT, only_mut
    
#This function makes the little gene diagrams for each transcript    
def gene_patches3(tx, tx_dict, ax):
    start, end, strand, CDS_start, CDS_end, exons = tx_info(tx, tx_dict)
    tx_patch = patches.Rectangle((start,0.8),end-start,0.02,edgecolor='0.1',facecolor='0.1')
    ax.add_patch(tx_patch)
    if exons is not None:
        exon_patches = []
        for exon_start, exon_stop in exons:
            exon_patches.append(patches.Rectangle((exon_start, 0.76),exon_stop-exon_start,0.09,edgecolor='0.1',facecolor='0.1'))
        for patch in exon_patches:
            ax.add_patch(patch)
    else:
        CDS_patch = patches.Rectangle((CDS_start, 0.76),CDS_end-CDS_start, 0.09,edgecolor='0.1',facecolor='0.1')
        ax.add_patch(CDS_patch)
    ax.get_yaxis().set_ticks([])
    return strand    

#Plot normalized WT and mutant
def plot_peak_differences(set1, set2, only1, only2, tx_dict, save_dir=None, use_max_y=False, tx_list=None):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    if tx_list is None:
        tx_list = [tx for tx in only1.keys() if len(only1[tx][1][0]) > 0]
    for tx in tx_list:     
        s1 = set1[tx][0]
        s2 = set2[tx][0]
        
        #Plot WT on top with markers for peaks only in WT
        fig, ax = plt.subplots(3, figsize=(12,6), sharex=True)
        ax[0].bar(s1.index, s1)
        ax[0].plot(only1[tx][1][0], only1[tx][1][1], 'o', color='orangered')
        
        #Plot mutant on the bottom (inverted) with markers for peaks only in mutant
        ax[1].bar(s2.index, s2, color='teal', edgecolor='teal')
        ax[1].plot(only2[tx][1][0], only2[tx][1][1], 'o', color='darkorange')

        fig.subplots_adjust(hspace=0)
        
        #Determine y_max - use the highest peak in either if use_max_y is True, otherwise use highest peak in lower
        if use_max_y is False:
            max_y = min([max(set1[tx][0]),max(set2[tx][0])])+0.2*min([max(set1[tx][0]),max(set2[tx][0])])
        else:
            max_y = max([max(set1[tx][0]),max(set2[tx][0])])+0.2*max([max(set1[tx][0]),max(set2[tx][0])])
        ax[0].set_ylim(0,max_y)
        ax[1].set_ylim(0,max_y)
        
        #Add the diagram of the gene
        ax3 = plt.subplot
        strand = gene_patches3(tx, tx_dict, ax[2])

        #Add title and plot
        ax[0].set_title(tx+" ("+strand+" strand)")
        ax[1].invert_yaxis()
        plt.show()
        
        #Save in sav_dir
        if save_dir is not None:
            fig.savefig(save_dir+tx+'.pdf', format='pdf')
        plt.clf()
        
#Pipeline to read bam files, call peaks, compare reproducibility between replicates and then compare wt and mut. Outputs plots as pdf. If you don't want all the plots in this directory, change save_dir
def compare_NETseq_pipeline(WTbam1, WTbam2, MUTbam1, MUTbam2, save_dir='./', gff3=None, tx_list=None):
    if gff3 is not None:
        tx_dict = build_transcript_dict(gff3)
    else: tx_dict = build_transcript_dict()
    tx_dict.update(build_centromere_dict())
    
    #Generate information from BAM files
    WT1 = NETseq_peaks(WTbam1, tx_dict, use_robust=False)
    WT2 = NETseq_peaks(WTbam2, tx_dict, use_robust=False)
    MUT1 = NETseq_peaks(MUTbam1, tx_dict, use_robust=False)
    MUT2 = NETseq_peaks(MUTbam2, tx_dict, use_robust=False)
    
    #Find reproducible peaks and average replicates
    WT_set, rep1, rep2 = compare_NETseq_peak_counts(WT1, WT2, from_reps=True)
    MUT_set, rep1, rep2 = compare_NETseq_peak_counts(MUT1, MUT2, from_reps=True)
    
    #Compare mutant and wild type
    WT_MUT_both, only_WT, only_MUT = compare_NETseq_peak_counts(WT_set, MUT_set)
    
    #Gather most abundant transcripts
    if tx_list is None:
        tx_list = []
        count = 0
        for tx, info in WT_set.iteritems():
            if len(info[0]) > 0:
                if sum(info[0])/len(info[0]) > .07 and 'RNA' not in tx and tx in new_tx_dict:
                    count += 1
                    spliced_list.append(tx)
        print count
        
    plot_peak_differences(WT_set, MUT_set, only_WT, only_MUT, tx_dict, save_dir=save_dir, tx_list=tx_list)