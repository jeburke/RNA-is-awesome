__author__ = 'jordanburke'
'''Usage: python Merge_lists.py <as many lists with the same first column header as you want> <output prefix> '''

import sys
import pandas as pd

def build_tables(file):
    fin = open(file, "r")
    df = pd.read_table(fin)
    print df.columns
    df.fillna(0)
    fin.close()
    return df

n = 2
while n < len(sys.argv)-1:
    file = sys.argv[n]
    df_merge = build_tables(sys.argv[1])
    df = build_tables(file)
    df_merge = pd.merge(df_merge,df,on=df.columns[0],left_index=True,how = 'left')
    n += 1
#print df_merge

df_merge_sorted = df_merge.sort(columns=df.columns[0])

fout = open("{0}_compared.txt".format(sys.argv[-1].split(".")[0]), "w")
print fout
df_merge_sorted.to_csv(fout, sep='\t')
fout.close()

