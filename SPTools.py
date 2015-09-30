__author__ = 'jordanburke'

import pandas



#################################################################
## Build dataframes from tables written by analyzeSp.py. Reads ##
## file containing Transcript column and columns with counts   ##
## for each type of read (exon, intron or total)               ##
#################################################################

def build_tables(file):
    fin = open(file, "r")
    df = pandas.read_table(fin)
    df.columns = df.columns.str.rstrip("_al.sorted.bam")
    df.fillna(0)
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
    merged = pandas.merge(df1,df2,on="Transcrip",left_index=True,how='left',suffixes=('_exons','_total'))
    merged = pandas.merge(merged,df3, on="Transcrip", left_index=True, how = 'left')
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
        if name.strip() == "Transcrip":
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
## Normalize counts in B samples to counts in W samples. Also  ##
## divides by transcript length. c1-c4 are read counts for     ##
## internal control RNA.                                       ##
#################################################################

def normalize_to_mature(total_counts, lengths, config):
    df1 = total_counts
    df2 = lengths
    df3 = config
    merged = pandas.merge(df1,df2, on="Transcrip", left_index=True, how = 'left')
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
        if name.strip() == "Transcrip":
            print "Reading table"
        elif n < len(configNames) and name == configNames[n]+"-B":
            print "Reading sample " +name
            c1 = controlValues[n][0]
            print "Control value 1 = " +c1
            c2 = controlValues[n][1]
            print "Control value 2 = " +c2
            c3 = controlValues[n][2]
            print "Control value 3 = " +c3
            untagged = controlValues[n][3]
            print "Untagged sample is " +untagged
            merged[name+"_norm"] = pandas.Series(((merged[name])/(merged[name.strip("B")+"W"])/merged[name.strip("Length")]), index = merged.index)
            merged[name+"_norm_untagged"] = pandas.Series((merged[name+"_norm"]/merged[untagged]), index = merged.index)
            n += 1
    merged.fillna(0)
    return merged

##################################################################
## Filter transcripts based on an input table. In this case I'm ##
## giving it a list of RNAi targets. Can be anything, 1 column  ##
## format with CNAG IDs                                         ##
##################################################################

def filter_transcripts(mergedtable, list):
    geneID = []
    fin = open(list, "r")
    df = pandas.DataFrame(mergedtable)
    df = df.set_index('Transcrip')
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

