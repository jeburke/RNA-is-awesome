__author__ = 'jordanburke'

import pandas
import math
import matplotlib.pyplot as plt

#################################################################
## Build dataframes from tables written by analyzeSp.py. Reads ##
## file containing Transcript column and columns with counts   ##
## for each type of read (exon, intron or total)               ##
#################################################################

def build_tables(file, header='infer', skiprows=None):
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
    merged = pandas.merge(df1_indexed,df2_indexed, left_index=True, right_index=True, how = 'left')
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
    merged = merged.set_index(("Transcript","Transcript"))
    merged = merged.fillna(0)
    print len(merged)
    return merged



#################################################################
## Normalize counts in B samples to counts in W samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA. Arguments are output from build_tables##
#################################################################

def normalize_to_mature(total_counts, lengths, config, cutoff=0):
    df1 = total_counts
    df2 = lengths
    df3 = config
    df1_indexed = df1.set_index(df1[("Transcript","Transcript")])
    df2_indexed = df2.set_index(df2[("Transcript","Transcript")])
    merged = pandas.merge(df1_indexed,df2_indexed, left_index=True, right_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df3.columns:
        print name
        configNames.append(name)
        s = pandas.Series(df3[name])
        l = s.tolist()
        controlValues.append(l)
    print controlValues
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
    RNAi_df = df[(df.index.isin(geneID))]
    print len(RNAi_df)
    return RNAi_df

#################################################################
## Convert output from normalize_counts to lists for plotting  ##
#################################################################

def get_ratios(normalizedResults, samplename, count_type, log=False):
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
            if y != 0:
                log_list.append(math.log10(y))
            else:
                log_list.append("NaN")
        values = log_list
    return values

###################################################################
## Plot ratios (usually log) of replicates for normalized counts ##
###################################################################


def scatter_plot(xvalues1, yvalues1, xvalues2, yvalues2, plot_title='3p ends/Total SP', legend1='All', legend2='Filtered'):
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    ax1.scatter(xvalues1, yvalues1, c='royalblue', label=legend1, alpha = 0.5, edgecolor='darkslateblue')
    ax1.scatter(xvalues2, yvalues2, c='coral', alpha=0.5, label=legend2, edgecolor='coral')
    ax1.set_xlabel("Replicate 1 (log10)")
    ax1.set_ylabel("Replicate 2 (log10)")
    xlim = ax1.get_xlim()
    ax1.set_ylim(xlim)
    ax1.set_xlim(xlim)
    ylim = ax1.get_ylim()
    ax1.legend(loc=4)
    ax1.set_title(plot_title)
    ax1.plot(xlim, ylim, ls="--", c=".3", lw=1.5)
    plt.show()