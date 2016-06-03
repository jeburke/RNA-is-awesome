__author__ = 'jordanburke'

import pandas
import math
import matplotlib.pyplot as plt
import numpy as np
import random
from beeswarm import *

#################################################################
## Convert output from normalize_counts to lists for plotting  ##
#################################################################

def get_ratios(dataFrame, samplename, count_type, log=False, base=10, T0only=False):
    values = []
    print len(dataFrame.columns)
    for name1, name2 in dataFrame.columns:
        #print name1, name2
        if name2 == count_type:
            #print name2
            if name1 == samplename:
                #print name1
                #print dataFrame
                s = pandas.Series(dataFrame[(name1,name2)])
                if T0only == True:
                    transcript_list = []
                    for row in dataFrame.iterrows():
                        if isinstance(row[0], tuple):
                            if row[0][0].endswith("T0"):
                                transcript_list.append(row[0][0])
                        else:
                            if row[0].endswith("T0"):
                                transcript_list.append(row[0])
                    s = s[s.index.get_level_values(0).isin(transcript_list)]
                #print s
                values = s.tolist()
            else:
                continue
    if log==True:
        log_list = []
        for y in values:
            if float(y) == float('Inf'):
                log_list.append(np.NaN)
            elif y == 0:
                log_list.append(np.NaN)
            else:
                log_list.append(math.log(y, base))
        values = log_list
    return values

##########################################################################
## Count splice sites or transcripts with at least some number of reads ##
##########################################################################

def count_above_cutoff(dataframe, sample_tuple, cutoff=5):
    counter1 = 0
    counter2 = 0

    for name1, name2 in dataframe.columns:
        if name1 == sample_tuple[0] and name2 == sample_tuple[1]:
            print name1, name2
            s = pandas.Series(dataframe[(name1,name2)])
            #print s
            for row in s.iteritems():
                #print row
                #print row.index
                if isinstance(row[0], tuple):
                    if row[0][0].endswith("T0"):
                        counter2 +=1
                        if row[1] > cutoff:
                            counter1 +=1                                   
                else:
                    if row[0].endswith("T0"):
                        #print row [0]
                        #print row [1]
                        counter2 +=1
                        if row[1] > cutoff:
                            counter1 +=1
                            
    if sample_tuple[1].endswith("prime"):
        print "Number of "+sample_tuple[1]+" splice sites with at least "+str(cutoff)+" reads: "+str(counter1)
        print "Number of "+sample_tuple[1]+" splice sites: "+str(counter2)
        print str(float(counter1)/float(counter2)*100.)+"% of splice sites"
    elif sample_tuple[1] == "Total":
        print "Number of transcripts with at least "+str(cutoff)+" reads: "+str(counter1)
        print "Number of "+sample_tuple[1]+" transcripts: "+str(counter2)
        print str(float(counter1)/float(counter2)*100.)+"% of transcripts"

#################################################################################                    
## Function that determines if values in the 1st sample are higher or lower    ##
##than values in the 2nd sample (2 fold)                                       ##
#################################################################################

def compare_samples(df, SampleTuple1, SampleTuple2, fold_change=10):
    count_type1 = SampleTuple1[1]
    count_type2 = SampleTuple2[1]
    samplename1 = SampleTuple1[0]
    samplename2 = SampleTuple2[0]
    
    for name1, name2 in df.columns:
        if name2 == count_type1 and name1 == samplename1:
            s = pandas.Series(df[(name1,name2)])
            dict1 =  s.to_dict()
        elif name2 == count_type2 and name1 == samplename2:
            s = pandas.Series(df[(name1,name2)])
            dict2 = s.to_dict()
            
    low_dict = {}
    low_counter = 0
    high_dict = {}
    high_counter = 0
    for name, value in dict1.iteritems():
        cutoff_high = float(dict2[name])*fold_change
        cutoff_low = float(dict2[name])/fold_change
        if name in dict2:
            if float(dict1[name]) > cutoff_high and dict2[name]!=0 and dict1[name]!=0:
                high_dict[name] = []
                high_dict[name].append(dict1[name])
                high_dict[name].append(dict2[name])
                high_dict[name].append(dict1[name]/dict2[name])
                high_counter += 1
            elif float(dict1[name]) < cutoff_low and dict2[name]!=0 and dict1[name]!=0:
                low_dict[name] = []
                low_dict[name].append(dict1[name])
                low_dict[name].append(dict2[name])
                low_dict[name].append(dict1[name]/dict2[name])
                low_counter += 1
    
    with open("{0}_{1}_comparison.tsv".format(SampleTuple1[0],SampleTuple2[0]), "w") as fout:
        fout.write("Transcript\t Intron\t Reads in "+SampleTuple1[0]+"\t Reads in "+SampleTuple2[0]+"\t Ratio\n")
        fout.write("Higher in "+SampleTuple1[0]+" than "+SampleTuple2[0]+"\n")
        for name, value in high_dict.iteritems():
            line_list = [name[0], str(name[1]), str(value[0]), str(value[1]), str(value[2]), "\n"]
            line = "\t".join(line_list)
            fout.write(line)
        fout.write("Lower in "+SampleTuple1[0]+" than "+SampleTuple2[0]+"\n")
        for name, value in low_dict.iteritems():
            line_list = [name[0], str(name[1]), str(value[0]), str(value[1]), str(value[2]), "\n"]
            line = "\t".join(line_list)
            fout.write(line)
    print str(high_counter)+" introns higher in "+SampleTuple1[0]+" than "+SampleTuple2[0]
    print str(low_counter)+" introns lower in "+SampleTuple1[0]+" than "+SampleTuple2[0]
    return (high_dict, low_dict)

def compare_dicts(dict1, dict2):
    both_dict = {}
    for name, values in dict1.iteritems():
        if name in dict2:
            both_dict[name] = values
    print len(both_dict)
    return both_dict
            
###################################################################
## Plot ratios (usually log) of replicates for normalized counts ##
###################################################################


def scatter_plot(xvalues1, yvalues1, xvalues2, yvalues2, plot_title='3p ends/Total SP', legend1='All', legend2='Filtered', xlabel='Replicate 1 (log10)', ylabel='Replicate 2 (log10)'):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(xvalues1, yvalues1, c='royalblue', label=legend1, alpha = 0.5, edgecolor='darkslateblue')
    ax1.scatter(xvalues2, yvalues2, c='coral', alpha=0.5, label=legend2, edgecolor='coral')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    print np.nanmin(xvalues1)
    print np.nanmax(xvalues1)
    xmax = np.nanmax(xvalues1)+0.5
    xmin = np.nanmin(xvalues1)-0.5
    ax1.set_ylim([xmin,xmax])
    ax1.set_xlim([xmin,xmax])
    ymax = ax1.get_ylim()
    ax1.legend(loc=4)
    ax1.set_title(plot_title)
    ax1.plot([xmin, xmax], [xmin,xmax], ls="--", c=".3", lw=1.5)
    plt.show()

#######################################################################################################################
## Scatter plot of normalized counts (log10 is the default).                                                         ##
## SampleTuple1 and 2 are tuples that contain the sample name and read type - e.g. ("CJB66D","Normalized to mature") ##
## DataFrames are from normalize_B_to_mature or normalize_AtoB.                                                      ##
## base is the base for log scale (if you want it to be 2, set base=2)                                               ##
#######################################################################################################################

def scatter_plot2(SampleTuple1, SampleTuple2, DataFrame1, DataFrame2, plot_title='3p ends/Total SP', legend1='All', legend2='Filtered', xlabel='Replicate 1 (log10)', ylabel='Replicate 2 (log10)', base=10, plot_both=True, log_both=True, scaling="auto"):
    print DataFrame1.columns.lexsort_depth
    #print DataFrame1.index.lexsort_depth
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    if plot_both==True and log_both==True:
        xvalues1 = get_ratios(DataFrame1, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues1 = get_ratios(DataFrame1, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        xvalues2 = get_ratios(DataFrame2, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues2 = get_ratios(DataFrame2, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        ax1.scatter(xvalues1, yvalues1, c='royalblue', label=legend1, alpha = 0.5, edgecolor='darkslateblue')
        ax1.scatter(xvalues2, yvalues2, c='coral', alpha=0.5, label=legend2, edgecolor='coral')
        ax1.legend(loc=4)
    elif plot_both==1 and log_both==True:
        xvalues1 = get_ratios(DataFrame1, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues1 = get_ratios(DataFrame1, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        ax1.legend(loc=4)
    elif plot_both==2 and log_both==True:
        xvalues2 = get_ratios(DataFrame2, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
        yvalues2 = get_ratios(DataFrame2, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
        ax1.scatter(xvalues2, yvalues2, c='0.3', alpha=1, label=legend2, edgecolor='0.2')
        ax1.legend(loc=4)
    elif plot_both==True and log_both==False:
        xvalues1 = get_ratios(DataFrame1, SampleTuple1[0], SampleTuple1[1], log=False)
        yvalues1 = get_ratios(DataFrame1, SampleTuple2[0], SampleTuple2[1], log=False)
        xvalues2 = get_ratios(DataFrame2, SampleTuple1[0], SampleTuple1[1], log=False)
        yvalues2 = get_ratios(DataFrame2, SampleTuple2[0], SampleTuple2[1], log=False)
        ax1.scatter(xvalues1, yvalues1, c='royalblue', label=legend1, alpha = 0.5, edgecolor='darkslateblue')
        ax1.scatter(xvalues2, yvalues2, c='royalblue', alpha=0.5, label=legend2, edgecolor='darkslateblue')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #print np.nanmin(xvalues1)
    #print np.nanmax(xvalues1)
    
    ax1.set_title(plot_title)
    if scaling == "auto":
        xmax = np.nanmax(xvalues1)+0.5
        xmin = np.nanmin(xvalues1)-0.5
        ax1.set_ylim([xmin,xmax])
        ax1.set_xlim([xmin,xmax])
        ymax = ax1.get_ylim()
        ax1.plot([xmin, xmax], [xmin, xmax], ls="--", c=".3", lw=1.5)
    plt.show()
    return fig1

def scatter_plot3(DataFrame_x, SampleList_x, DataFrame_y, SampleList_y, plot_title='Plot', legend1='All', legend2='Filtered', xlabel='Replicate 1 (log10)', ylabel='Replicate 2 (log10)', base=10):
    a=0
    xvalues = []
    while a < len(SampleList_x):
        xvalues.append(get_ratios(DataFrame_x, SampleList_x[a][0], SampleList_x[a][1], log=True, base=base, T0only=True))
        a += 1
    b=0
    yvalues = []
    while b < len(SampleList_y):
        yvalues.append(get_ratios(DataFrame_y, SampleList_y[b][0], SampleList_y[b][1], log=True, base=base, T0only=True))
        b += 1
    print len(xvalues)
    print len(yvalues)
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    n = 0
    color_list = ["royalblue","coral","lightsteelblue"]
    while n < len(yvalues):
        if len(xvalues) == len(yvalues):   
            print len(xvalues[n])
            print len(yvalues[n])
            ax1.scatter(xvalues[n], yvalues[n], c=color_list[n], label=legend1, alpha = 0.5, edgecolor=color_list[n])
            n += 1
        elif len(xvalues) == 1:
            print len(xvalues[0])
            print len(yvalues[n])
            ax1.scatter(xvalues[0], yvalues[n], c=color_list[n], label=legend1, alpha = 0.5, edgecolor=color_list[n])
            n += 1
        else:
            print "X and Y variables do not match"
            break
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #print np.nanmin(xvalues[1])
    #print np.nanmax(xvalues[1])
    #xmax = np.nanmax(xvalues[1])+0.5
    #xmin = np.nanmin(xvalues[1])-0.5
    #ax1.set_ylim([xmin,xmax])
    #ax1.set_xlim([xmin,xmax])
    #ymax = ax1.get_ylim()
    ax1.legend(loc=4)
    ax1.set_title(plot_title)
    ax1.plot([xmin, xmax], [xmin,xmax], ls="--", c=".3", lw=1.5)
    plt.show()
    return fig1


#######################################################################################################################
## Bar chart of averaged replicates of normalized counts (error bars are the range).                                 ##
## df is DataFrame from normalize_B_to_mature or normalize_AtoB                                                      ##
## sample1 and 2 are sample names (e.g. CJB66D). These should be in level0 of your DataFrame Columns                 ##
## CNAG_lists are the lists of transcripts you want to plot. 1 and 2 will be different colors. You only need to give ##
## CNAG_list1                                                                                                        ##
#######################################################################################################################

def bar_chart(df, read_type, sample1, sample2, CNAG_list1, CNAG_list2=None, plot_title="Stuff", ylabel="Stuff"):
    if CNAG_list2 is not None:
        CNAG_list = CNAG_list1 + CNAG_list2
    else:
        CNAG_list = CNAG_list1
    gene_list = []
    values1 = []
    values2 = []
    #Check if the dataframe index has only transcript names or both transcripts and exons
    if type(df.index[0]) == str:
        for name1, name2 in df:
            if name1 == sample1 and name2 == read_type:
                for gene in CNAG_list:
                    if gene in df.index:
                        values1.append(df[(name1,name2)][gene])
                        gene_list.append(gene)
            elif name1 == sample2 and name2 == read_type:
                for gene in CNAG_list:
                    if gene in df.index:
                        values2.append(df[(name1,name2)][gene])
    else:
        for gene in CNAG_list:
            for transcript in df.iterrows():
                if transcript[0][0] == gene:
                    values1.append(df[(sample1,read_type)][(transcript[0][0], transcript[0][1])])
                    values2.append(df[(sample2,read_type)][(transcript[0][0], transcript[0][1])])
                    gene_list.append(str(transcript[0][0])+"-"+str(transcript[0][1]))
    n=0
    avg_values = []
    errors = []
    while n < len(values1):
        avg_values.append((values1[n]+values2[n])/2)
        errors.append(abs(values1[n]-values2[n])/2)
        n += 1
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    n_groups = len(avg_values)
    index = np.arange(n_groups)
    bar_width = 0.5
    error_config = {'ecolor': '0.3'}
    ax1.bar(index, avg_values, bar_width,alpha=0.9,color='darkslateblue',yerr=errors,error_kw=error_config)
    plt.xticks(index, gene_list, rotation='vertical')
    plt.xlabel('Transcript')
    plt.ylabel(ylabel)
    plt.title(plot_title)
    print ax1.get_children()
    for child in ax1.get_children()[4:4+len(CNAG_list1)]:
        child.set_color('coral')
    plt.show()
    return fig1


#############################################################################
## Get a random set of transcripts from a file with a list of transcripts. ##
## Pretty self explanatory                                                 ##
#############################################################################

def random_transcripts_from_list(file_with_list, number_of_transcripts):
    gene_list = []
    fin = open(file_with_list, "r")
    for line in fin:
        if line.startswith("CNAG"):
            gene = line.split("\t")[0].strip()
            if gene[-2:-1] == "T":
                gene_list.append(gene)
            else:
                gene_list.append(gene+"T0")
    random_list = random.sample(gene_list, number_of_transcripts)
    return random_list

def beeswarm_plot(DataFrame_list, SampleTuple_list, DataFrame2=None, base=10, color_list="blue", select_random=False):
    print SampleTuple_list
    print len(SampleTuple_list)
    values_list = []
    name_list = []
    median_list = []
    a=0
    while a < len(DataFrame_list):
        n=0
        while n < len(SampleTuple_list):
            print SampleTuple_list[n]
            values = get_ratios(DataFrame_list[a], SampleTuple_list[n][0], SampleTuple_list[n][1], log=True, base=base)
            median_list.append(np.median(np.array(values)))
            if select_random != False:
                values = random.sample(values, select_random)
            values_list.append(values)
            name_list.append(SampleTuple_list[n][0])
            n += 1
        a += 1
    print median_list
    bs, ax = beeswarm(values_list, method="swarm", labels=name_list, col=color_list)