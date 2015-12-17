__author__ = 'jordanburke'
'''Usage: python Merge_lists.py <as many lists with the same first column header as you want> <output prefix> '''

import sys
import pandas as pd

def build_tables(file):
    fin = open(file, "r")
    transcripts = []
    info = []
    header = []
    for line in fin:
        line_list = (line.strip()).split("\t")
        if line.startswith("CNAG"):
            transcripts.append(line_list[0]+"-"+line_list[1])
            a = 2
            info_list = []
            while a < len(line_list):
                info_list.append(line_list[a])
                a+=1
            info.append(info_list)
        else:
            b = 1
            while b < len(line_list):
                print len(line_list)
                print b
                header.append(line_list[b])
                b += 1
    print len(transcripts)
    print len(info)
    print header
    ds_info = pd.Series(info)
    transcript_dict = dict(zip(transcripts,ds_info))
#    print transcript_dict['CNAG_01026T0-8']
    df = pd.DataFrame()
    df = df.from_dict(transcript_dict, orient='index')
    df.columns = header[0:len(df.columns)]
    print df
    return df

n = 2
df2 = build_tables(sys.argv[1])
while n < len(sys.argv)-1:
    file = sys.argv[n]
    df1 = build_tables(file)
    print df1.index.name
    df2 = pd.merge(df2, df1, right_index=True, left_index=True)
    print df2.index.name
    n += 1

df2_sorted = df2.sort(columns=df2.columns[0])

fout = open("{0}_compared.txt".format(sys.argv[-1].split(".")[0]), "w")
print fout
df2_sorted.to_csv(fout, sep='\t')
fout.close()


