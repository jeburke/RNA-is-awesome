__author__ = 'jordanburke'

import pandas
import math
import matplotlib.pyplot as plt
import numpy as np
import random

#################################################################
## Build dataframes from tables written by analyzeSp.py. Reads ##
## file containing Transcript column and columns with counts   ##
## for each type of read (exon, intron or total)               ##
#################################################################

def build_tables(file, header='infer', skiprows=None, low_memory=False):
    #For junction reads and cleavages use header=[0,1] and skiprows=[2]. For all other tables, use default values.
    fin = open(file, "r")
    df = pandas.read_table(fin, header=header, skiprows=skiprows)
    columns_list = []
    for column in df.columns:
        if header != 'infer':
            if column[0] == "Unnamed: 0_level_0":
                columns_list.append(("Transcript","Transcript"))
            elif column[0] == "Unnamed: 1_level_0":
                columns_list.append(("Exon","Exon"))
            else:
                columns_list.append((column[0].split("_")[0], column[1]))
        else:
            if column == "Unnamed: 0":
                columns_list.append(("Transcript","Transcript"))
            elif column == "Transcript":
                columns_list.append(("Transcript","Transcript"))
            else:
                columns_list.append((column.split("_")[0],"Total"))
    df.columns = pandas.MultiIndex.from_tuples(columns_list)
    df = df.fillna(0)
    return df

##################################################################
## Build dataframes from tables written by analyzeSp.py.        ##
## This function concatenates the Transcript and Intron columns ##
## with a "-" to compare to other lists, such as splice site    ##
## scores.                                                      ##
##################################################################

def build_intron_tables(file):
    fin = open(file, "r")
    transcripts = []
    info = []
    header = []
    for line in fin:
        line_list = (line.strip()).split("\t")
        if line.startswith("CNAG"):
            if line_list[0][-2:] == "T0":
                transcripts.append(line_list[0]+"-"+line_list[1])
                a = 2
                info_list = []
                while a < len(line_list):
                    info_list.append(line_list[a])
                    a+=1
                info.append(info_list)
            else:
                transcripts.append(line_list[0]+"T0-"+str(int(line_list[1])+1))
                a = 2
                info_list = []
                while a < len(line_list):
                    info_list.append(line_list[a])
                    a+=1
                info.append(info_list)
        elif line.startswith("Transcript"):
            b = 1
            while b < len(line_list):
                header.append(line_list[b].strip())
                b += 1
    ds_info = pandas.Series(info)
    transcript_dict = dict(zip(transcripts,ds_info))
    df = pandas.DataFrame()
    df = df.from_dict(transcript_dict, orient='index')
    df.columns = header[1:len(df.columns)+1]
    return df

#################################################################
## Normalize counts in A samples to counts in B samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA. Arguments are output from build_tables##
#################################################################

def normalize_AtoB(feature_counts, total_counts, lengths, config, cutoff=0):
    df1 = feature_counts
    df2 = total_counts
    df3 = lengths
    df4 = config
    df1_indexed = df1.set_index(df1[("Transcript","Transcript")])
    df2_indexed = df2.set_index(df2[("Transcript","Transcript")])
    df3_indexed = df3.set_index(df3[("Transcript","Transcript")])
    merged = pandas.merge(df1_indexed, df2_indexed, left_index=True, right_index=True, how = 'left')
    merged = pandas.merge(merged, df3_indexed, left_index=True, right_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df4.columns:
        configNames.append(name)
        s = pandas.Series(df4[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    for name, read_type in merged.columns:
        names.append(name)
        if name.strip() == "Transcript":
            print "Reading table"
        elif n < len(configNames) and name == configNames[n][0]+"-A":
            print "Reading sample " +name
            c1 = controlValues[n][0]
            print "Control value 1 = " +str(c1)
            c2 = controlValues[n][1]
            print "Control value 2 = " +str(c2)
            c3 = controlValues[n][2]
            print "Control value 3 = " +str(c3)
            merged[(name.split("-")[0],"5prime Normalized")] = pandas.Series((((merged[(name,"5prime")])/c1)/((merged[(name.strip("-A")+"-B","Total")])/(merged[(name.strip("Length"),"Total")]*c2))), index = merged.index)
            merged[(name.split("-")[0],"3prime Normalized")] = pandas.Series((((merged[(name,"3prime")])/c1)/((merged[(name.strip("-A")+"-B","Total")])/(merged[(name.strip("Length"),"Total")]*c2))), index = merged.index)
            if cutoff > 0:
                merged = merged[merged[(name.split("-")[0]+"-B","Total")] > cutoff]
            n += 1
    if merged.columns[1] == ("Exon","Exon"):
        merged.set_index([merged.columns[0], merged.columns[1]], inplace=True)
    else:
        merged = merged.set_index(("Transcript","Transcript"))
    merged = merged.fillna(0)
    print "5prime Normalized; 3prime Normalized"
    print len(merged)
    return merged



#################################################################
## Normalize counts in B samples to counts in W samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA. Arguments are output from build_tables##
#################################################################

def normalize_B_to_mature(total_counts, lengths, config, cutoff=0):
    df1 = total_counts
    df2 = lengths
    df3 = config
    df1_indexed = df1.set_index(df1[("Transcript","Transcript")])
    df2_indexed = df2.set_index(df2[("Transcript","Transcript")])
    merged = pandas.merge(df1_indexed,df2_indexed, left_index=True, right_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df3.columns:
        configNames.append(name)
        s = pandas.Series(df3[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    for name, count_type in merged.columns:
        names.append(name)
        for config_name, sample_type in df3.columns:
            if name == config_name+"-B":
                print "Reading sample " +name
                c1 = df3[(config_name, sample_type)][0]
                print "Control value 1 = " +str(c1)
                c2 = df3[(config_name, sample_type)][1]
                print "Control value 2 = " +str(c2)
                c3 = df3[(config_name, sample_type)][2]
                print "Control value 3 = " +str(c3)
                merged[(name.split("-")[0],"Normalized to mature")] = pandas.Series(((merged[(name,"Total")]/c2)/(merged[(name.split("-")[0]+"-W","Total")]/c3)), index = merged.index)
                if cutoff > 0:
                    merged = merged[merged[(name.split("-")[0]+"-B","Total")] > cutoff]
    merged = merged.fillna(0)
    print "Normalized to mature"
    print len(merged)
    return merged


####################################################################
## Filter transcripts based on an input table. In this case I'm   ##
## giving it a list of RNAi targets. Can be anything, 1 column    ##
## format with CNAG IDs. mergedtable is from one of the normalize ##
## functions. list_file is the name of a files containing the list##
## of CNAGs you want to filter by.                                ##
####################################################################

def filter_transcripts_by_cnag(mergedtable, list_file):
    geneID = []
    fin = open(list_file, "r")
    df = pandas.DataFrame(mergedtable)
    for line in fin:
        if line.startswith("CNAG"):
            gene = line.strip()
            geneID.append(gene)
            geneID.append(gene+"T0")
            geneID.append(gene+"T1")
            geneID.append(gene+"T2")
            geneID.append(gene+"T3")
            geneID.append(gene+"T4")
            geneID.append(gene+"T5")
    RNAi_df = pandas.DataFrame(columns=df.columns)
    print len(df.index)
    if type(df.index[0]) == str:
        RNAi_df = df[(df.index.isin(geneID))]
    else:
        index_list = []
        for transcript in df.iterrows():
            if transcript[0][0] in geneID:
                RNAi_df = RNAi_df.append(df.loc[transcript[0]])
                index_list.append(transcript[0])
                RNAi_df.index = pandas.MultiIndex.from_tuples(index_list)
    print len(RNAi_df)
    return RNAi_df

#################################################################
## Convert output from normalize_counts to lists for plotting  ##
#################################################################

def get_ratios(normalizedResults, samplename, count_type, log=False, base=10):
    values = []
    for name1, name2 in normalizedResults.columns:
        if name2 == count_type:
            if name1 == samplename:
                s = pandas.Series(normalizedResults[(name1, name2)])
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

def scatter_plot2(SampleTuple1, SampleTuple2, DataFrame1, DataFrame2, plot_title='3p ends/Total SP', legend1='All', legend2='Filtered', xlabel='Replicate 1 (log10)', ylabel='Replicate 2 (log10)', base=10):
    xvalues1 = get_ratios(DataFrame1, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
    print SampleTuple1[0]
    print SampleTuple1[1]
    yvalues1 = get_ratios(DataFrame1, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
    xvalues2 = get_ratios(DataFrame2, SampleTuple1[0], SampleTuple1[1], log=True, base=base)
    yvalues2 = get_ratios(DataFrame2, SampleTuple2[0], SampleTuple2[1], log=True, base=base)
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
    for child in ax1.get_children()[1:1+len(CNAG_list1)]:
        child.set_color('coral')
    plt.show()


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

    
