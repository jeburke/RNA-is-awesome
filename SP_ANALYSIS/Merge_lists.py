__author__ = 'jordanburke'

import sys
import pandas as pd

def build_tables(file):
    fin = open(file, "r")
    df = pd.read_table(fin)
    df.columns = df.columns.str.strip()
    df.fillna(0)
    fin.close()
    return df


n = 2
while n < len(sys.argv):
    df_merge = build_tables(sys.argv[1])
    "df"+n = build_tables(sys.argv[n])
    df_merge = pd.merge(df_merge,"df"+n, on=[0], left_index=True, how = 'left')
    n += 1

fout = open("{0}_compared.txt".format(sys.argv[-1].split(".")[0]), "w")
fout.write(pandas.DataFrame.to_csv(df_merge, sep='\t'))
fout.close()

