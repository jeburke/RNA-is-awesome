import json
import pysam
import pandas as pd
import numpy as np
import sys
import os

def open_gff3(gff3_file):
    gff3 = pd.read_csv(gff3_file, sep='\t', comment='#', 
            names=['chromosome','source','type','start','end','x','strand','y','ID'])
    gff3.loc[:,'name'] = gff3['ID'].str.split(';').str[0].str.lstrip('ID=')
    return gff3

def polA_ds(read, chrom, fa, frag_size):
    polA = False
    if not read.is_reverse:
        ds_seq = fa[chrom][read.reference_end:read.reference_end+frag_size]
        if 'AAAAAAAA' in ds_seq: polA = True
    else:
        ds_seq = fa[chrom][read.reference_start-frag_size:read.reference_start]
        if 'TTTTTTTT' in ds_seq: polA = True
    return polA

def open_fa(fa):
    if fa.endswith('json'):
        with open(fa) as f:
            fa_dict = json.load(f)
    else:
        with open(fa) as f:
            fa_dict = {}
            for line in f:
                if line.startswith('>'):
                    contig = line.lstrip('>').split()[0]
                    fa_dict[contig] = ''
                else:
                    fa_dict[contig] = fa_dict[contig]+line.strip()
    return fa_dict

def count_reads(bam_file, gff3_file, extend=0, report_unique=False, detect_polA=False, fa=None, frag_size=100):
    gff3 = open_gff3(gff3_file)    
    bam = pysam.Samfile(bam_file)
    if detect_polA:
        fa = open_fa(fa)

    counts = pd.DataFrame(index=gff3[gff3['type'] == 'gene']['name'], columns = ['All'])
    counted_reads = {}
    polA_counter = 0
    ext_counter = 0
    print "Counting all reads in each feature..."
    for ix, gene in gff3[gff3['type'] == 'gene'].iterrows():
        counted_reads[gene['name']] = set()
        
        if gene['strand'] == '+': 
            start = gene['start']
            end = gene['end'] + extend
        elif gene['strand'] == '-': 
            start = gene['start'] - extend
            if start < 0: start = 0
            end = gene['end']

        region = bam.fetch(gene['chromosome'], start, end)
        
        for read in region:
            if detect_polA:
                if polA_ds(read, gene['chromosome'], fa, frag_size):
                    polA_counter += 1
                    continue
                    
            if read.reference_start in range(gene['start'],gene['end']+1) and read.reference_end in range(gene['start'],gene['end']+1):
                if read.is_reverse and gene['strand'] == '+':
                    counted_reads[gene['name']].add(read)
                if not read.is_reverse and gene['strand'] == '-':
                    counted_reads[gene['name']].add(read)
            else:
                if gene['strand'] == '+' and read.reference_start - extend <= gene['end'] and not read.is_reverse:
                    counted_reads[gene['name']].add(read)
                    ext_counter += 1
                elif gene['strand'] == '-' and read.reference_end + extend >= gene['start'] and read.is_reverse:
                    counted_reads[gene['name']].add(read)
                    ext_counter += 1
        counts.loc[gene['name'],'All'] = len(counted_reads[gene['name']])
        
    if extend > 0:
        print str(ext_counter)+' reads detected by extending gene'
        
    if detect_polA:
        print '\n'+str(polA_counter)+' reads neighboring polyA sequence removed'
    
    if report_unique:
        red_counter = 0
        print "\nDetecting non-unique reads..."
        for gene in counts.index:
            overlap = set()
            for other_gene in list(set(counts.index).difference(gene)):
                overlap.union(counted_reads[gene].intersection(counted_reads[other_gene]))
            counts.loc[gene,'Unique'] = counts.loc[gene,'All']-len(overlap)
            red_counter += len(overlap)
        print str(red_counter)+' non-unique reads found'
            
    return counts

def main():
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print '''Usage:
       Required arguments
       ------------------
       --gff3 : str
                   gff3 file location
       
       **Provide either --directory or --bam_file**
       --directory : str
                   Directory containing multiple bam files (will count reads in all bam files with indexes)
       --bam_file : str
                   Location of single bam file to process
       
       Optional arguments
       ------------------
       --report-unique : bool, default `False`
                   In addition to counting all reads, also report reads that are uniquely assigned to one feature
       --extend-reads : int, default 0
                   Extend the 5p end of the read n bp - helpful when annotation does not include 3' UTRs
       --remove_pA_adj_reads : bool, defuault `False`
                   Remove reads that neighbor polyA stretches (control for false priming)
                   If True, also provide --fasta and --insert_size
       --fasta : str
                   Fasta file for genome, provide if remove_pA_adj_reads is True
       --insert_size : int, default 100
                   Approximate size of library inserts - used to determine how far away to look for polA tracts
                        
    '''
        return None
    
    gff3_file = sys.argv[sys.argv.index('--gff3')+1]
    
    if '--directory' in sys.argv:
        directory = sys.argv[sys.argv.index('--directory')+1]
        if not directory.endswith('/'): directory = directory+'/'
        bam_files = []
        for file in os.listdir(directory):
            if file.endswith('.bam') and file+'.bai' in os.listdir(directory):
                bam_files.append(directory+file)
    elif '--bam_file' in sys.argv:
        bam_file = sys.argv[sys.argv.index('--bam_file')+1]
        bam_files = [bam_file]
    else:
        print "Must provide directory or bam_file"
        return None
    
    report_unique = False
    if '--report-unique' in sys.argv:
        report_unique = True
        
    extend = 0
    if '--extend-reads' in sys.argv:
        extend = int(sys.argv[sys.argv.index('--extend-reads')+1])
        
    detect_polA = False
    fa = None
    frag_size = 100
    if '--remove_pA_adj_reads' in sys.argv:
        detect_polA = True
        try:
            fasta = sys.argv[sys.argv.index('--fasta')+1]
        except ValueError:
            print "Must provide fasta file to identify reads neighboring pA"
            return None
        
        try:
            frag_size = int(sys.argv[sys.argv.index('--insert_size')+1])
        except ValueError:
            frag_size = 100
    
    final_df = None
    for bam_file in bam_files:
        df = count_reads(bam_file, gff3_file, extend=extend, report_unique=report_unique, detect_polA=detect_polA, fa=fasta, frag_size=frag_size)
        bam_name = bam_file.split('/')[-1].split('.bam')[0]
        df.columns = pd.MultiIndex.from_product([[bam_name], df.columns])
        if final_df is None:
            final_df = df
        else:
            final_df = final_df.merge(df, right_index=True, left_index=True, how='outer')
    
    name = 'QuantSeq_counts.csv'
    if '--name' in sys.argv:
        name = sys.argv[sys.argv.index('--name')+1]+'.csv'
    
    final_df.to_csv(name)
        
if __name__ == "__main__":
    main()