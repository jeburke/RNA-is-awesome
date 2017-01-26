__author__ = 'jordanburke'

import pandas
import math
import matplotlib.pyplot as plt
import numpy as np
import random
from beeswarm import *
from scipy.stats import ks_2samp
import re
import operator
import seaborn as sns
sns.set_style('white')

#################################################################
## Convert output from normalize_counts to lists for plotting  ##
#################################################################

def get_ratios(dataFrame, samplename, count_type, log=False, base=10, T0only=False):
    values = []
    for name1, name2 in dataFrame.columns:
        if name2 == count_type:
            if name1 == samplename:
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
        if plot_both == True or plot_both == 1:
            xmax = np.nanmax(xvalues1)+0.5
            xmin = np.nanmin(xvalues1)-0.5
        elif plot_both == 2:
            xmax = np.nanmax(xvalues2)+0.5
            xmin = np.nanmin(xvalues2)-0.5
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
    
def histogram(df1, SampleTuple1, df2=None, SampleTuple2=None, xlabel="Intermediate levels", ylabel="Number of introns", bins=100, legend1='All', legend2=None):
    x1 = get_ratios(df1, SampleTuple1[0], SampleTuple1[1], log=True, base=10)
    x1 = [x for x in x1 if str(x) != 'nan']
    fig1 = plt.figure()
    ax = fig1.add_subplot(111)
    ax.hist(x1, bins=bins, color='royalblue', edgecolor='royalblue', alpha=0.8, label=legend1 )
    if df2 is not None and SampleTuple2 is None:
        x2 = get_ratios(df2, SampleTuple1[0], SampleTuple1[1], log=True, base=10)
        x2 = [x for x in x2 if str(x) != 'nan']
        ax.hist(x2, bins=bins, color='coral', edgecolor='coral', alpha=0.8, label=legend2)
    elif SampleTuple2 is not None and df2 is None:
        y1 = get_ratios(df1, SampleTuple2[0], SampleTuple2[1], log=True, base=10)
        y1 = [x for x in y1 if str(x) != 'nan']
        ax.hist(y1, bins=bins, color='coral', edgecolor='coral', alpha=0.8, label=legend2)
    elif df2 is not None and SampleTuple2 is not None:
        y2 = get_ratios(df2, SampleTuple2[0], SampleTuple2[1], log=True, base=10)
        y2 = [x for x in y2 if str(x) != 'nan']
        ax.hist(y2, bins=bins, color='coral', edgecolor='coral', alpha=0.8, label=legend2)
    if df2 is not None or SampleTuple2 is not None:    
        ax.legend(loc=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    return fig1
    
##Input is the outuput from the Score_splice_sites.py script (Splice_site_strengths.tsv), can also take score_dict directly
def bin_transcripts(input_file, bin_size, bin_by="length", df_list=None):
    score_dict = {}
    if type(input_file) == str:
        with open(input_file, "r") as fin:
            for line in fin:
                if line.startswith("Transcript"): 
                    continue
                else:
                    columns = re.split(r'\t', line)
                    transcript = columns[0]
                    intron = columns[1]
                    five_p_score = float(columns[2])
                    three_p_score = float(columns[3])
                    length = int(columns[4])
                    key = (transcript, int(intron)+1)
                    if bin_by == "length" and length > 30:
                        score_dict[key] = length
                    elif bin_by == "length" and length <= 30:
                        continue
                    elif bin_by == "five_p_score":
                        score_dict[key] = five_p_score
                    elif bin_by == "three_p_score":
                        score_dict[key] = three_p_score
                    else: print "I don't recognize the bin_by value"
    elif type(input_file) == dict:
        score_dict = input_file
    if df_list is not None:
        score_dict = {key: score_dict[key] for key in score_dict if key in df_list }
    sorted_dict = sorted(score_dict.items(), key=operator.itemgetter(1))

    n = 0
    bin1 = []
    values1 = []
    bin2 = []
    values2 = []
    all_scores = []
    for transcript, score in sorted_dict:
        all_scores.append(score)
        if n < bin_size:
            bin1.append(transcript)
            values1.append(score)
            n += 1
        elif n >= bin_size and n < len(sorted_dict)-bin_size:
            n += 1
        else:
            bin2.append(transcript)
            values2.append(score)
            n += 1
    print "Mean intron length or splice site score:"
    print "%0.1f" % np.mean(all_scores)
    print "Mean intron length or splice site score for bottom bin:"
    print "%0.1f" % np.mean(values1)
    print "Mean intron length or splice site score for top bin:"
    print "%0.1f" % np.mean(values2)
    return (bin1, bin2)
        
    
def bin_by_position(df, first_last=True, first_other=False, other_last=False):
    df_list = df.index.tolist()
    intron_dict = {}
    for intron in df_list:
        if intron[0] not in intron_dict:
            intron_dict[intron[0]]=[]
        intron_dict[intron[0]].append(intron[1])
    
    bin1 = []
    bin2 = []
    for transcript, introns in intron_dict.iteritems():
        introns = sorted(introns)
        first = introns[0]
        last = introns[-1]
        if first_last == True:
            bin1.append((transcript, first))
            bin2.append((transcript, last))
        elif first_other == True:
            bin1.append((transcript, first))
            i=1
            for i in range(len(introns)):
                bin2.append((transcript, introns[i]))
        elif other_last == True:
            bin2.append((transcript, last))
            i=0
            for i in range(len(introns)-1):
                bin1.append((transcript, introns[i]))
    if len(bin1) != len(bin2):
        if len(bin1) < len(bin2):
            bin2 = random.sample(bin2, len(bin1))
        elif len(bin1) > len(bin2):
            bin1 = random.sample(bin1, len(bin2))
    return bin1, bin2
    
def cumulative_function(df, score_file, bin_size, SampleTuple1, bin_by=False, bin_tuple=None, SampleTuple2=None, xlabel="Intermediate levels", ylabel="Number of introns", title="CDF", plot_type="CDF", bins=100):
    #df = df[df[SampleTuple1] > 0]
    if bin_by == 'position':
        transcript_list_1, transcript_list_2 = bin_by_position(df, first_last=False, first_other=False, other_last=True)
    elif bin_by != False:
        df_list = df.index.tolist()
        transcript_list_1, transcript_list_2 = bin_transcripts(score_file, bin_size, bin_by, df_list=df_list)
    elif bin_tuple is not None:
        transcript_list_1= bin_tuple[0]
        transcript_list_2= bin_tuple[1]
    new_df1 = pandas.DataFrame(columns=df.columns)
    new_df2 = pandas.DataFrame(columns=df.columns)

    new_df1 = df[df.index.map(lambda x: x in transcript_list_1)]
    new_df2 = df[df.index.map(lambda x: x in transcript_list_2)]
    
    fig1 = plt.figure(figsize=(8, 6), dpi=600)
    ax1 = fig1.add_subplot(111)
    if plot_type == "CDF":
        x1 = get_ratios(new_df1, SampleTuple1[0], SampleTuple1[1], log=False)
        x2 = get_ratios(new_df2, SampleTuple1[0], SampleTuple1[1], log=False)
        values1, base1 = np.histogram(x1, bins=bins)
        values2, base2 = np.histogram(x2, bins=bins)
        cumulative1 = np.cumsum(values1)
        base1 = np.insert(base1, 0, 0)
        cumulative1 = np.insert(cumulative1, 0, 0)
        ax1.plot(base1[:-1], cumulative1, c='blue', linewidth=3.0, label="Low 1")
        cumulative2 = np.cumsum(values2)
        base2 = np.insert(base2, 0, 0)
        cumulative2 = np.insert(cumulative2, 0, 0)
        ax1.plot(base2[:-1], cumulative2, c='orangered', linewidth=3.0, label="High 1")
        if SampleTuple2 is not None:
            y1 = get_ratios(new_df1, SampleTuple2[0], SampleTuple2[1], log=False)
            values3, base3 = np.histogram(y1, bins=1000)
            y2 = get_ratios(new_df2, SampleTuple2[0], SampleTuple2[1], log=False)
            values4, base4 = np.histogram(y2, bins=1000)
            cumulative3 = np.cumsum(values3)
            base3 = np.insert(base3, 0, 0)
            cumulative3 = np.insert(cumulative3, 0, 0)
            ax1.plot(base3[:-1], cumulative3, color='lightskyblue', linewidth=3.0, label="Low 2")
            cumulative4 = np.cumsum(values4)
            base4 = np.insert(base4, 0, 0)
            cumulative4 = np.insert(cumulative4, 0, 0)
            ax1.plot(base4[:-1], cumulative4, color='coral', linewidth=3.0, label="High 2")
            ax1.legend(loc=4)
    
    elif plot_type == "PDF":
        x1 = get_ratios(new_df1, SampleTuple1[0], SampleTuple1[1], log=True)
        x2 = get_ratios(new_df2, SampleTuple1[0], SampleTuple1[1], log=True)
        x1 = [x for x in x1 if str(x) != 'nan']
        x2 = [x for x in x2 if str(x) != 'nan']
        #x1 = [1e-15 if str(x) == 'nan' else x for x in x1]
        #x2 = [1e-15 if str(x) == 'nan' else x for x in x2]
        print "Zero values removed:"
        print str(bin_size-len(x2))+" from high bin"
        print str(bin_size-len(x1))+" from low bin"
        #ax1.hist(x2, color='coral', edgecolor='coral', label="High 1", bins=bins, alpha=0.9)
        #ax1.hist(x1, color='royalblue', edgecolor='royalblue', label="Low 1", bins=bins, alpha=0.5)
        sns.distplot(x2, color='orangered', label="High", ax=ax1, bins=bins)
        sns.distplot(x1, color='royalblue', label="Low", ax=ax1, bins=bins)
        ax1.legend(loc=1)
        if SampleTuple2 is not None:
            y1 = get_ratios(new_df1, SampleTuple2[0], SampleTuple2[1], log=True)
            y2 = get_ratios(new_df2, SampleTuple2[0], SampleTuple2[1], log=True)
            y1 = [x for x in y1 if str(x) != 'nan']
            y2 = [x for x in y2 if str(x) != 'nan']
            ax1.hist(y2, color='orangered',  label="High 2", bins=bins)
            ax1.hist(y1, color='lightskyblue', label="Low 2", bins=bins)

    print "KS_statistic, p_value for replicate 1: "
    ks_x1 = get_ratios(new_df1, SampleTuple1[0], SampleTuple1[1], log=False)
    ks_x2 = get_ratios(new_df2, SampleTuple1[0], SampleTuple1[1], log=False)
    for value in ks_2samp(ks_x1, ks_x2):
        print "%0.1e" % value
    if SampleTuple2 is not None:
        print "KS_statistic, p_value for replicate 2: "
        ks_y1 = get_ratios(new_df1, SampleTuple2[0], SampleTuple2[1], log=False)
        ks_y2 = get_ratios(new_df2, SampleTuple2[0], SampleTuple2[1], log=False)
        for value in ks_2samp(ks_y1, ks_y2):
            print "%0.1e" % value
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    #ax1.set_ylim([0.9*len(new_df1),len(new_df1)+len(new_df1)*0.01])
    #ax1.set_ylim([0*len(new_df1),len(new_df1)+len(new_df1)*0.01])
    xmax = np.nanmax([np.nanmax(x1),np.nanmax(x2)])
    ax1.set_xlim([-1,xmax+.05*xmax])
    #ax1.set_xlim([-2, 50])
    ax1.set_title(title)
    plt.show()
    return fig1


##################################################################################
## Functions for normalizing data and creating a CDF that actually adds up to 1 ##
## Takes a list of data                                                         ##
##################################################################################

def normalize(data):
    data = [float(x)/sum(data) for x in data]
    return data

def cdf_values(data, bins='auto'):
    if bins=='auto':
        bins = len(data)/10
    values, base = np.histogram(data, bins=bins)
    values = normalize(values)
    cumulative = np.cumsum(values)
    base = np.insert(base, 0, min(data))
    cumulative = np.insert(cumulative, 0, 0)
    return cumulative, base

def cdf_for_n_lists(list_of_lists, label_list=None, color_list=None, x_title='Lengths', y_title='Fraction of introns'):
    fig = plt.figure(figsize=(8, 6), dpi=600)
    ax = fig.add_subplot(111)
    if label_list is None:
        label_list = range(len(list_of_lists))
    if color_list is None:
        color_list = ['0.3','cornflowerblue','orangered','0.7','limegreen','mediumvioletred']
    
    n = 0 
    for n in range(len(list_of_lists)):
        cumulative, base = cdf_values(list_of_lists[n])
        ax.plot(base[:-1], cumulative, c=color_list[n], linewidth=3.0, label=label_list[n])
        
    #ax.legend(bbox_to_anchor=(1, 0), loc='lower left', fontsize=14)
    ax.legend(fontsize=14)
    plt.ylabel(y_title, fontsize=16)
    plt.xlabel(x_title, fontsize=16)
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)
    plt.show()
    
    return fig
        
def sns_distplot_n(list_of_lists, label_list=None, color_list=None, x_title='Lengths', y_title='Franction of introns'):
    fig = plt.figure(figsize=(8, 6), dpi=600)
    ax = fig.add_subplot(111)
    if label_list is None:
        label_list = range(len(list_of_lists))
    if color_list is None:
        color_list = ['0.3','cornflowerblue','orangered','0.7','limegreen','mediumvioletred']
    
    n = 0 
    for n in range(len(list_of_lists)):
        ax = sns.distplot(list_of_lists[n])
        
    