__author__ = 'jordanburke'

import pandas
import math
import matplotlib.pyplot as plt

#################################################################
## Build dataframes from tables written by analyzeSp.py. Reads ##
## file containing Transcript column and columns with counts   ##
## for each type of read (exon, intron or total)               ##
#################################################################

def build_tables(file):
    fin = open(file, "r")
    df = pandas.read_table(fin)
    new_columns = []
    for column in df.columns:
        new_columns.append(column.split("_")[0])
    df.columns = new_columns
    #df.columns = df.columns.str.rstrip("_sorted.bam")
    df.fillna(0)
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
## internal control RNA.                                       ##
#################################################################

def normalize_counts(exon_counts,total_counts,lengths,config):
    df1 = exon_counts
    df2 = total_counts
    df3 = lengths
    df4 = config
    merged = pandas.merge(df1,df2,on="Transcript",left_index=True,how='left',suffixes=('_exons','_total'))
    merged = pandas.merge(merged,df3, on="Transcript", left_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df4.columns:
        configNames.append(name)
        s = pandas.Series(df4[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    for name in merged.columns:
        names.append(name)
        if name.strip() == "Transcript":
            print "Reading table"
        elif n < len(configNames) and name == configNames[n]+"-A_exons":
            c1 = controlValues[n][0]
            print c1
            c2 = controlValues[n][1]
            print c2
            merged[name+"_norm"] = pandas.Series((merged[name]/c1)/(merged[name.strip("A_exons")+"B_total"]/c2), index = merged.index)
            n += 1
    merged.fillna(0)
    return merged

#################################################################
## Normalize counts in A samples to counts in B samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA.                                       ##
#################################################################

def normalize_AtoB(feature_counts, total_counts, lengths, config):
    df1 = feature_counts
    df2 = total_counts
    df3 = lengths
    df4 = config
    merged = pandas.merge(df1,df2, on="Transcript", left_index=True, how = 'left', suffixes = ('_end', '_total'))
    merged = pandas.merge(merged, df3, on="Transcript", left_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df4.columns:
        configNames.append(name)
        s = pandas.Series(df4[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    for name in merged.columns:
        names.append(name)
        if name.strip() == "Transcript":
            print "Reading table"
        elif n < len(configNames) and name == configNames[n]+"-A_end":
            print "Reading sample " +name
            c1 = controlValues[n][0]
            print "Control value 1 = " +str(c1)
            c2 = controlValues[n][1]
            print "Control value 2 = " +str(c2)
            c3 = controlValues[n][2]
            print "Control value 3 = " +str(c3)
            merged[name+"_norm"] = pandas.Series((((merged[name])/c1)/((merged[name.strip("-A_end")+"-B_total"])/(merged[name.strip("Length")]*c2))), index = merged.index)

            n += 1
    merged = merged.set_index("Transcript")
    merged = merged.fillna(0)
    return merged


def normalize_AtoB_filter_by_counts(feature_counts, total_counts, lengths, config, cutoff):
    df1 = feature_counts
    df2 = total_counts
    df3 = lengths
    df4 = config
    merged = pandas.merge(df1,df2, on="Transcript", left_index=True, how = 'left', suffixes = ('_end', '_total'))
    merged = pandas.merge(merged, df3, on="Transcript", left_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df4.columns:
        configNames.append(name)
        s = pandas.Series(df4[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    for name in merged.columns:
        names.append(name)
        if name.strip() == "Transcript":
            print "Reading table"
        elif n < len(configNames) and name == configNames[n]+"-A_end":
            print "Reading sample " +name.split("_")[0]
            c1 = controlValues[n][0]
            print "Control value 1 = " +str(c1)
            c2 = controlValues[n][1]
            print "Control value 2 = " +str(c2)
            c3 = controlValues[n][2]
            print "Control value 3 = " +str(c3)
            merged[name+"_norm"] = pandas.Series((((merged[name])/c1)/((merged[name.strip("-A_end")+"-B_total"])/(merged[name.strip("Length")]*c2))), index = merged.index)
            merged_filtered = merged[merged[name.split("-")[0]+"-B_total"] > cutoff]
            n += 1
    merged_filtered = merged_filtered.set_index("Transcript")
    merged_filtered = merged_filtered.fillna(0)
    return merged_filtered


#################################################################
## Normalize counts in B samples to counts in W samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA.                                       ##
#################################################################

def normalize_to_mature(total_counts, lengths, config):
    df1 = total_counts
    df2 = lengths
    df3 = config
    merged = pandas.merge(df1,df2, on="Transcript", left_index=True, how = 'left')
    configNames = []
    controlValues = []
    for name in df3.columns:
        configNames.append(name)
        s = pandas.Series(df3[name])
        l = s.tolist()
        controlValues.append(l)
    names = []
    n = 0
    for name in merged.columns:
        names.append(name)
        if name.strip() == "Transcript":
            print "Reading table"
        elif n < len(configNames) and name == configNames[n]+"-B":
            print "Reading sample " +name
            c1 = controlValues[n][0]
            print "Control value 1 = " +str(c1)
            c2 = controlValues[n][1]
            print "Control value 2 = " +str(c2)
            c3 = controlValues[n][2]
            print "Control value 3 = " +str(c3)
            merged[name+"_norm"] = pandas.Series(((merged[name])/(merged[name.strip("B")+"W"])/merged[name.strip("Length")]), index = merged.index)
            n += 1
    merged = merged.fillna(0)
    return merged


##################################################################
## Filter transcripts based on an input table. In this case I'm ##
## giving it a list of RNAi targets. Can be anything, 1 column  ##
## format with CNAG IDs                                         ##
##################################################################

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
    return RNAi_df

#################################################################
## Convert output from normalize_counts to lists for plotting  ##
#################################################################

def get_ratios(normalizedResults, a):
    values = []
    n = 0
    for name in normalizedResults.columns:
        if name[-5:] == "_norm":
            if n == a:
                s = pandas.Series(normalizedResults[name])
                values = s.tolist()
                n += 1
            else:
                n += 1
    return values

#################################################################
## Take log10 of data to make plots prettier                   ##
#################################################################

def log_ratios(ratio_list):
    log_list = []
    for y in ratio_list:
        if float(y) != 0:
            y = float(y)
            log_list.append(math.log10(y))
        elif float(y) == 0:
            log_list.append("NaN")
    return log_list

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