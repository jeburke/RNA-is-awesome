"""
Usage:
    Compare_filter_CGmap_both.py  (--output_file_path=PATH)(--minimum_methylation=PERCENT)(--min_counts1=READ_NUMBER)(--min_counts2=READ_NUMBER)(--input1_file_path=path_to_CGmap_file)(--input2_file_path=path_to_CGmap_file)
"""
import pandas as pd
import sys
from docopt import docopt
import numpy as np
import matplotlib
import os

def compare_CGmap(CG_file1, CG_file2, out_path):
    print "Comparing CGmap files..."
    CG1 = pd.read_csv(CG_file1, sep='\t', header=None,
                     names=['chromosome','base','position','W-strand pair','C-strand pair','fraction methylated','counts_meth','counts'])
    CG2 = pd.read_csv(CG_file2, sep='\t', header=None,
                     names=['chromosome','base','position','W-strand pair','C-strand pair','fraction methylated','counts_meth','counts'])
    
    both = CG1.merge(CG2, right_on=['chromosome','base','position','W-strand pair','C-strand pair'], 
                     left_on=['chromosome','base','position','W-strand pair','C-strand pair'], how='outer',suffixes=(' 1',' 2'))
    
    both.index = both['chromosome'].str.cat(both['position'].apply(int).apply(str), sep=':')
    both = both.dropna(how='any')

    return both

def filter_CGmap(df, min_counts1, min_counts2, min_meth, output_file_path, CG_file1, CG_file2):
    print "Positions in unfiltered: "+str(len(df))
    print len(df[df['counts 1'] >= min_counts1])
    print len(df[df['counts 2'] >= min_counts2])
    
    # Initial filter for number of counts and fraction methylated to minimize data frame size while iterating
    df = df[(df['counts 1'] >= min_counts1) & (df['counts 2'] >= min_counts2)]
    df = df[(df['fraction methylated 1'] >= min_meth-0.1) | (df['fraction methylated 2'] >= min_meth-0.1)]

    print "After initial filter: "+str(len(df))
    df_array = df.as_matrix()
    # columns are [chromosome, base, position, W-strand pair, C-strand pair,
    # fraction methylated 1, counts_meth 1, counts 1, fraction methylated 2,
    # counts_meth 2, counts 2]
    
    col_ix = {'chromosome':0,'base':1,'position':2,'W-strand pair':3,'C-strand pair':4,
       'fraction methylated 1':5,'counts_meth 1':6,'counts 1':7,
       'fraction methylated 2':8,'counts_meth 2':9,'counts 2':10}
    
    new_ix = []
    col_dict = {}
    for col in df.columns:
        col_dict[col] = []
    col_dict['Filter 1'] = []
    col_dict['Filter 2'] = []
    
    # Iterate through matrix to check that each site has a sequential partner
    n=0
    nonseq = 0
    while n in range(len(df_array)-1):
        if n%1000==0:
            print n
        r = df_array[n]
        next_r = df_array[n+1]
        if next_r[0]==r[0] and next_r[2]-r[2] == 1:
            # Check if either data set falls below the cutoffs
            flag1=True
            flag2=True
            # Check if first sample has requisite minimum counts and methylation
            if (r[7]<min_counts1 or next_r[7]<min_counts1 or r[5]<min_meth or 
                next_r[5]<min_meth or str(r[7]) == 'nan'):
                flag1=False
                
            # Check if second sample has requisite minimum counts and methylation    
            if (r[10]<min_counts2 or next_r[10]<min_counts2 or r[8]<min_meth or 
                next_r[8]<min_meth or str(r[10]) == 'nan'):
                flag2=False

            new_ix.append(df.index[n])
            new_ix.append(df.index[n+1])
            
            # Build new columns from data in original matrix according to new index
            for key in col_dict:
                if 'Filter' not in key:
                    col_dict[key].append(r[col_ix[key]])
                    col_dict[key].append(next_r[col_ix[key]])
            
            # Categorize pairs of CGs as pass or fail based on filters above
            if flag1:
                col_dict['Filter 1'].append('Pass')
                col_dict['Filter 1'].append('Pass')
            else:
                col_dict['Filter 1'].append('Fail')
                col_dict['Filter 1'].append('Fail')
            if flag2:
                col_dict['Filter 2'].append('Pass')
                col_dict['Filter 2'].append('Pass')
            else:
                col_dict['Filter 2'].append('Fail')
                col_dict['Filter 2'].append('Fail')
            n += 2
        else:
            nonseq += 1
            n += 1
    
    new_df = pd.DataFrame(columns=df.columns, index=new_ix)
    print "Number of rows without sequential next row: "+str(nonseq)
    
    for col, data in col_dict.iteritems():
        new_df.loc[:,col] = data

    new_df = new_df[(new_df['Filter 1'] == 'Pass') & (new_df['Filter 2'] == 'Pass')]
    print "Sites in both CGmaps: "+str(len(new_df))
    
    new_df['fraction methylated'] = new_df[['fraction methylated 1','fraction methylated 2']].mean(axis=1)
    new_df['counts_meth'] = new_df[['counts_meth 1','counts_meth 2']].min(axis=1)
    new_df['counts'] = new_df[['counts 1','counts 2']].min(axis=1)
    
    out_cols = ['chromosome','base','position','W-strand pair','C-strand pair','fraction methylated','counts_meth','counts']
    CG_name = '{0}/{1}_{2}_both.CGmap'.format(output_file_path,CG_file1.split('/')[-1].split('.CGmap')[0],
                                                 CG_file2.split('/')[-1].split('.CGmap')[0])
    print CG_name
    new_df[out_cols].to_csv(CG_name, header=False, index=False, sep='\t')
    
    return new_df

def make_bedgraphs(df, CG_file1, CG_file2, output_file_path, out_path):
    
    bg1 = df[['chromosome','position','fraction methylated']]
    bg1.loc[:,'position'] = bg1['position'].apply(int)
    bg1.loc[:,'position-1'] = bg1['position']-1
    bg1 = bg1[['chromosome','position-1','position','fraction methylated']]
    bg1_name = '{0}/{1}_{2}_both.bedgraph'.format(output_file_path,CG_file1.split('/')[-1].split('.CGmap')[0],
                                                 CG_file2.split('/')[-1].split('.CGmap')[0])
    print bg1_name
    bg1.to_csv(bg1_name, sep='\t', header=False, index=False)

def main():
    arguments=docopt(__doc__,version="Compare_filter_CGmap 0.1")
    for k,v in arguments.iteritems():
        print k+': '+v
    input1_file_path=arguments["--input1_file_path"]
    input2_file_path=arguments["--input2_file_path"]
    output_file_path=arguments["--output_file_path"]
    min_meth=float(arguments["--minimum_methylation"])/100
    min_counts1=int(arguments["--min_counts1"])
    min_counts2=int(arguments["--min_counts2"])
    
    out_file = output_file_path.split('/')[-1]
    out_path = output_file_path.split('/')[0]
    both = compare_CGmap(input1_file_path, input2_file_path, out_path)
    filt = filter_CGmap(both, min_counts1, min_counts2, min_meth, output_file_path,input1_file_path, input2_file_path)
    make_bedgraphs(filt, input1_file_path, input2_file_path, output_file_path, out_path)

if __name__ == "__main__":
    main()


