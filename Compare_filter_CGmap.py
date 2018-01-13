"""
Usage:
    Compare_filter_CGmap.py  (--output_file_path=PATH)(--minimum_methylation=PERCENT)(--min_counts1=READ_NUMBER)(--min_counts2=READ_NUMBER)(--input1_file_path=path_to_CGmap_file)(--input2_file_path=path_to_CGmap_file)
"""
import pandas as pd
import sys
from docopt import docopt
import numpy as np
import matplotlib
import os

def compare_CGmap(CG_file1, CG_file2):
    CG1 = pd.read_csv(CG_file1, sep='\t', header=None,
                     names=['chromosome','base','position','W-strand pair','C-strand pair','fraction methylated','counts_meth','counts'])
    CG2 = pd.read_csv(CG_file2, sep='\t', header=None,
                     names=['chromosome','base','position','W-strand pair','C-strand pair','fraction methylated','counts_meth','counts'])
    
    both = CG1.merge(CG2, right_on=['chromosome','base','position','W-strand pair','C-strand pair'], 
                     left_on=['chromosome','base','position','W-strand pair','C-strand pair'], how='outer',suffixes=(' 1',' 2'))
    
    both.index = both['chromosome'].str.cat(both['position'].apply(int).apply(str), sep=':')
    one_not_two = both[both['fraction methylated 2'].apply(str) == 'nan'].index
    two_not_one = both[both['fraction methylated 1'].apply(str) == 'nan'].index
    
    print "\nSites in 1 but not in 2:"
    print len(one_not_two)
    
    print "Sites in 2 but not in 1:"
    print len(two_not_one)
    
    CG1.index = CG1['chromosome'].str.cat(CG1['position'].apply(int).apply(str), sep=':')
    CG1[CG1.index.isin(one_not_two)].to_csv(CG_file1.split('.CGmap')[0]+'_difference.CGmap', header=False, index=False, sep='\t')
    CG2.index = CG2['chromosome'].str.cat(CG2['position'].apply(int).apply(str), sep=':')
    CG2[CG2.index.isin(two_not_one)].to_csv(CG_file2.split('.CGmap')[0]+'_difference.CGmap', header=False, index=False, sep='\t')
    
    bg1 = CG1[CG1.index.isin(one_not_two)][['chromosome','position','fraction methylated']]
    bg1['position-1'] = bg1['position']-1
    bg1 = bg1[['chromosome','position-1','position','fraction methylated']]
    bg1.to_csv(CG_file1.split('.CGmap')[0]+'_difference.bedgraph', sep='\t', header=False, index=False)
    
    bg2 = CG2[CG2.index.isin(two_not_one)][['chromosome','position','fraction methylated']]
    bg2['position-1'] = bg2['position']-1
    bg2 = bg2[['chromosome','position-1','position','fraction methylated']]
    bg2.to_csv(CG_file2.split('.CGmap')[0]+'_difference.bedgraph', sep='\t', header=False, index=False)
    
    return both

def filter_CGmap(df, min_counts1, min_counts2, min_meth, output_file_path):
    print "Peaks in unfiltered: "+str(len(df))
    new_ix = []
    col_dict = {}
    for col in df.columns:
        col_dict[col] = []
    col_dict['Filter 1'] = []
    col_dict['Filter 2'] = []
    
    n=0
    nonseq = 0
    while n in range(len(df)):
        r = df.iloc[n]
        r_name = r.name
        next_r = df.iloc[n+1]
        next_r_name = next_r.name
        if next_r['chromosome']==r['chromosome'] and next_r['position']-r['position'] == 1:
            # Check if either data set falls below the cutoffs
            flag1=True
            flag2=True
            r = r.append(pd.Series(['Pass','Pass'], index=['Filter 1','Filter 2']))
            next_r = next_r.append(pd.Series(['Pass','Pass'], index=['Filter 1','Filter 2']))
            
            if (r['counts 1']<min_counts1 or next_r['counts 1']<min_counts1 or r['fraction methylated 1']<min_meth or 
                next_r['fraction methylated 1']<min_meth or str(r['counts 1']) == 'nan'):
                flag1=False
                r.loc['Filter 1'] = 'Fail'
                next_r.loc['Filter 1'] = 'Fail'
            if (r['counts 2']<min_counts2 or next_r['counts 2']<min_counts2 or r['fraction methylated 2']<min_meth or 
                next_r['fraction methylated 2']<min_meth or str(r['counts 2']) == 'nan'):
                flag2=False
                r.loc['Filter 2'] = 'Fail'
                next_r.loc['Filter 2'] = 'Fail'
            
            # If either passes the filter
            if flag1 or flag2:
                new_ix.append(r_name)
                new_ix.append(next_r_name)
                for key in col_dict:
                    col_dict[key].append(r[key])
                    col_dict[key].append(next_r[key])

            n += 2
            #print len(new_ix)
        else:
            #print "No sequential position: "+df.iloc[n].name
            nonseq += 1
            n += 1
    
    new_df = pd.DataFrame(columns=df.columns, index=new_ix)
    print "Peaks in filtered: "+str(len(new_df))
    print "Number of rows without sequential next row: "+str(nonseq)
    
    for col, data in col_dict.iteritems():
        new_df.loc[:,col] = data
    
    new_df['ratio'] = new_df['fraction methylated 2']/new_df['fraction methylated 1']
    new_df['log2 ratio'] = new_df['ratio'].apply(np.log2)
    new_df.loc[new_df[new_df['counts 1'].apply(str) == 'nan'].index,'log2 ratio'] = 7.0
    new_df.loc[new_df[new_df['counts 2'].apply(str) == 'nan'].index,'log2 ratio'] = -7.0
    
    new_df.to_csv(output_file_path)
    
    return new_df

def make_bedgraphs(df, CG_file1, CG_file2, output_file_path):
    # A - one_not_two and two_not_one after filter
    one_not_two = df[df['fraction methylated 2'].apply(str) == 'nan']
    two_not_one = df[df['fraction methylated 1'].apply(str) == 'nan']
    
    bg1 = one_not_two[['chromosome','position','fraction methylated 1']]
    bg1.loc[:,'position'] = bg1['position'].apply(int)
    bg1.loc[:,'position-1'] = bg1['position']-1
    bg1 = bg1[['chromosome','position-1','position','fraction methylated 1']]
    bg1.to_csv(CG_file1.split('.CGmap')[0]+'_filt_difference.bedgraph', sep='\t', header=False, index=False)
    
    bg2 = two_not_one[['chromosome','position','fraction methylated 2']]
    bg2.loc[:,'position'] = bg2['position'].apply(int)
    bg2.loc[:,'position-1'] = bg2['position']-1
    bg2 = bg2[['chromosome','position-1','position','fraction methylated 2']]
    bg2.to_csv(CG_file2.split('.CGmap')[0]+'_filt_difference.bedgraph', sep='\t', header=False, index=False)
    
    # B - log ratio bedgraph with all non-NaN points
    bg3 = df[df['log2 ratio'].apply(str) != 'nan'][['chromosome','position','log2 ratio']]
    bg3.loc[:,'position'] = bg3['position'].apply(int)
    bg3.loc[:,'position-1'] = bg3['position']-1
    bg3 = bg3[['chromosome','position-1','position','log2 ratio']]
    bg3.to_csv(output_file_path.split('.CGmap')[0].split('csv')[0]+'_all_log2ratio.bedgraph', sep='\t', header=False, index=False)
    
    # C - log ratio bedgraph with only passing in both data sets
    bg4 = df[(df['log2 ratio'].apply(str) != 'nan') & (df['Filter 1'] == 'Pass') & (df['Filter 2'] == 'Pass')]
    bg4.loc[:,'position'] = bg4['position'].apply(int)
    bg4.loc[:,'position-1'] = bg4['position']-1
    bg4 = bg4[['chromosome','position-1','position','log2 ratio']]
    bg4.to_csv(output_file_path.split('.CGmap')[0].split('.csv')[0]+'_filt_log2ratio.bedgraph', sep='\t', header=False, index=False)

def main():
    arguments=docopt(__doc__,version="filter_CGmap 0.1")
    for k,v in arguments.iteritems():
        print k+': '+v
    input1_file_path=arguments["--input1_file_path"]
    input2_file_path=arguments["--input2_file_path"]
    output_file_path=arguments["--output_file_path"]
    min_meth=float(arguments["--minimum_methylation"])/100
    min_counts1=int(arguments["--min_counts1"])
    min_counts2=int(arguments["--min_counts2"])
    
    both = compare_CGmap(input1_file_path, input2_file_path)
    filt = filter_CGmap(both, min_counts1, min_counts2, min_meth, output_file_path)
    make_bedgraphs(filt, input1_file_path, input2_file_path, output_file_path)

if __name__ == "__main__":
    main()


