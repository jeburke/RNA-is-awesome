'''Script for calling 5' and 3' splice sites and relative strengths of each sites 

Usage: python Call_splice_site_motifs.py

This will generate an output spreadsheet with relative splice-site strength for each transcript in the genome'''

import sys
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneUtility
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/')
import SPTools as SP
import collections

def read_juncbase_output(juncbase_output):
    if type(juncbase_output) == list:
        df_list = []
        for output in juncbase_output:
            with open(juncbase_output, 'r') as fin:
                df = pd.read_csv(fin, sep='\t')
                df_list.append(df)
        junc_df = pd.concat(df_list)
        junc_df.replace(to_replace='NA', value=np.NaN, inplace=True)
    elif type(juncbase_output) == str:
        with open(juncbase_output, 'r') as fin:
            junc_df = pd.read_csv(fin, sep='\t')
    
    sample_list = []
    junc_df['contained in'] = None
    n_max = len(junc_df.columns)
    n = 11
    while n < n_max:
        if junc_df.columns[n].startswith('avg'): continue
        elif junc_df.columns[n] == 'contained in': n += 1
        elif n == 11:
            sample_list.append(junc_df.columns[n])
            sample_list.append(junc_df.columns[n+1])
            junc_df['avg_'+junc_df.columns[n]+'+'+junc_df.columns[n+1]] = junc_df[junc_df.columns[n]]
            wt = junc_df.columns[n]
            n += 2
        else:
            junc_df['avg_'+junc_df.columns[n]+'+'+junc_df.columns[n+1]] = junc_df[[(junc_df.columns[n]),(junc_df.columns[n+1])]].mean(axis=1, skipna=False)
            sample_list.append(junc_df.columns[n])
            sample_list.append(junc_df.columns[n+1])
            if str(junc_df.columns[n]) == 'nan' or str(junc_df.columns[n+1]) == 'nan':
                print junc_df['avg_'+junc_df.columns[n]+'+'+junc_df.columns[n+1]]
            if n != 11:
                junc_df[junc_df.columns[n]+'+'+junc_df.columns[n+1]] = ''
            n += 2
    
    print junc_df.columns
    filt_df = pd.DataFrame(columns=junc_df.columns, index = junc_df.index)
    n = 0
    for n in range(len(junc_df)):
        sample_str = ''
        flag = False
        row = junc_df.iloc[[n]]
        a = 0
        for column in row.columns:
            if column.startswith('avg'):
                if a == 0:
                    wt_name = column
                elif a > 0:
                    if str(junc_df[column][n]) == 'nan':
                        junc_df[column.split('avg_')[1]] = np.NaN
                    else:
                        flag = True
                        sample_str = sample_str+column.split('avg_')[1]+', '
                        if junc_df[column][n] > junc_df[wt_name][n]:
                            junc_df[column.split('avg_')[1]] = 'Up'
                        elif junc_df[column][n] < junc_df[wt_name][n]:
                            junc_df[column.split('avg_')[1]] = 'Down'
                a += 1                   
#        if flag == True:
        filt_df.iloc[[n]] = junc_df.iloc[[n]]
        filt_df.loc[n,'contained in'] = sample_str

    filt_df.replace(to_replace=np.NaN, value='NA', inplace=True)
    filt_df = filt_df[filt_df['strand'] != 'NA']
    filt_df = filt_df.reset_index()
    
    n = 0
    for n in range(len(filt_df)):
        for column in filt_df.columns:
            if column in sample_list:
                if filt_df[column][n] == 'NA':
                    filt_df.loc[n, column] = filt_df[wt][n]
                    
    
    filt_df['coord_1'] = filt_df['exclusion_junctions'].str.split(':').str.get(1).str.split('-').str.get(0)
    filt_df['coord_2'] = filt_df['exclusion_junctions'].str.split(':').str.get(1).str.split('-').str.get(1).str.split(';').str.get(0)
    
    print sample_list
    return (filt_df, sample_list)

def make_fasta_dict(fasta_file):
    fasta_dict = {}
    letters = ['A','C','G','T','N']
    n = 0
    with open(fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                n += 1
                line_list = line.split(' ')
                chr_num = line_list[0].strip()[1:]
                print chr_num
                fasta_dict[chr_num] = str()
            elif line[0] in letters and "chr"+str(n) in fasta_dict:
                fasta_dict["chr"+str(n)] = fasta_dict["chr"+str(n)]+line.strip()

    fasta_dict = collections.OrderedDict(sorted(fasta_dict.items()))
    return fasta_dict


##Input: dictionary of coordinates by transcript (transcript:[list of coordinates][list of coordinates])
def get_sequence(coord_dict, gff3_file, fasta_file):
    transcript_dict = SP.build_transcript_dict(gff3_file)
    if type(fasta_file) is str:
        fasta_dict = make_fasta_dict(fasta_file)
    else:
        fasta_dict = fasta_file
    
    seq_dict = {}
    counter5 = 0
    counter3 = 0
    other = 0
    for transcript, coord_sets in coord_dict.iteritems():
        seq_dict[transcript] = []
        chrom = transcript_dict[transcript][3]
        strand = transcript_dict[transcript][2]
        for coord in coord_sets[0]:
            seq_type = 'other'
            if strand == "+":
                sequence = fasta_dict[chrom][(coord-9):(coord+11)]
            elif strand == "-":
                sequence = fasta_dict[chrom][(coord-10):(coord+10)]
                sequence = SP.reverse_complement(sequence)

            if sequence[10:12] == 'GT' or sequence[10:12] == 'GC': 
                seq_type = "5'"
                counter5 += 1
            seq_dict[transcript].append((sequence, seq_type))
     
        for coord in coord_sets[1]:
            seq_type = 'other'
            if strand == "+":
                sequence = fasta_dict[chrom][(coord-9):(coord+11)]
            elif strand == "-":
                sequence = fasta_dict[chrom][(coord-10):(coord+10)]
                sequence = SP.reverse_complement(sequence)
                
            if sequence[8:10] == 'AG': 
                seq_type = "3'"
                counter3 += 1
            seq_dict[transcript].append((sequence, seq_type))
    
    print str(counter5)+" 5' splice sites"
    print str(counter3)+" 3' splice sites"
    
    return seq_dict
        
def get_junction_sequence(df, gff3_file, fasta_file):
    df = df.sort_values('chr', axis=0)
    
    #transcript_dict[transcript] = [start, end, strand, chromosome, CDS start, CDS end]
    transcript_dict = SP.build_transcript_dict(gff3_file)

    #splice_dict[transcipt] = [[5'sites][3'sites]]
    splice_dict, flag = SP.list_splice_sites(gff3_file)
    
    #fasta_dict[chr] = sequence
    if type(fasta_file) is str:
        fasta_dict = make_fasta_dict(fasta_file)
    else:
        fasta_dict = fasta_file

    transcript_by_chr = {}
    for transcript, coords in transcript_dict.iteritems():
        chromosome = coords[3]
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(transcript)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(transcript)

    df['Gene'] = "Unknown"
    df['intron'] = "Middle"
    df['sequence1'] = ''
    df['sequence2'] = ''
    df['intron sequence'] = 'No sequence here'

    n = 0
    for n in range(len(df)):
        coord1 = int(df['coord_1'][n].strip())
        coord2 = int(df['coord_2'][n].strip())
        chrom = df['chr'][n].strip()
        strand = df['strand'][n].strip()
        transcripts = transcript_by_chr[chrom]

        for transcript in transcripts:
            tx_strand = transcript_dict[transcript][2]
            start = transcript_dict[transcript][0]
            stop = transcript_dict[transcript][1]
            
            if strand == tx_strand and coord1 >= start and coord2 <= stop:
                df.loc[n,'Gene'] = transcript
               
        if strand == '+':
            sequence1 = fasta_dict[chrom][(coord1-11):(coord1+9)]
            sequence2 = fasta_dict[chrom][(coord2-10):(coord2+10)]
            all_seq = fasta_dict[chrom][(coord1-1):coord2]
        elif strand == '-':
            sequence1 = fasta_dict[chrom][(coord2-10):(coord2+10)]
            sequence1 = SP.reverse_complement(sequence1)
            sequence2 = fasta_dict[chrom][(coord1-11):(coord1+9)]
            sequence2 = SP.reverse_complement(sequence2)
            all_seq = fasta_dict[chrom][(coord1-1):coord2]
            all_seq = SP.reverse_complement(all_seq)
        
        df.loc[n,'sequence1'] = sequence1
        df.loc[n,'sequence2'] = sequence2
        df.loc[n,'intron sequence'] = all_seq

    for transcript in transcripts:
        if transcript in df['Gene'].tolist():
            tx_df = df[df['Gene'] == transcript]
            s = tx_df['coord_1']
            min_idx = s.idxmin()
            first = int(s.min())
            print transcript_dict[transcript][2]
            print first
            max_idx = s.idxmax()
            last = int(s.max())
            print last
        
            if first == last:
                df.loc[min_idx,'intron'] = 'Only'
            else:
                if transcript_dict[transcript][2] == '+':
                    df.loc[min_idx,'intron'] = 'First'
                    df.loc[max_idx,'intron'] = 'Last'
                elif transcript_dict[transcript][2] == '-':
                    df.loc[min_idx,'intron'] = 'Last'
                    df.loc[max_idx,'intron'] = 'First'
            
            for index, coord_1 in s.iteritems():
                if df['intron'][index] == 'Middle':
                    if coord_1 in range(first-10, first+10):
                        df_idx = s[s == coord_1].index[0]
                        if transcript_dict[transcript][2] == '+':
                            df.loc[df_idx, 'intron'] = 'First'
                        elif transcript_dict[transcript][2] == '-':
                            df.loc[df_idx, 'intron'] = 'Last'
                    elif coord_1 in range(last-10, last+10):
                        df_idx = s[s == coord_1].index[0]
                        if transcript_dict[transcript][2] == '+':
                            df.loc[df_idx, 'intron'] = 'Last'
                        elif transcript_dict[transcript][2] == '-':
                            df.loc[df_idx, 'intron'] = 'First'
                
    df = df[df['contained in'] != '']
    df = df.reset_index()
    return df

def mirrored(listtree):
    return [(mirrored(x) if isinstance(x, list) else x) for x in reversed(listtree)]

def generate_consensus_matrix(fasta_dict):
    #Populate gene dictionary and build genome
    gene_list_as_dict = GeneUtility.GeneDict()
    gene_list_as_dict.PopulateFromFile_new()
    genome = fasta_dict
    print genome.keys()

    #First generate a consensus matrix for 5' and 3' splice site, where 1st row is A counts, second row is C, third row is T, fourth row is G.
    pos_matrix_5prime = np.zeros([4,20])
    pos_matrix_3prime = np.zeros([4,20])

    counter1 = 0
    counter2 = 0

    for cnag in gene_list_as_dict:
        counter2 += 1
        gene = gene_list_as_dict.FindGeneByCNAG(cnag)
        introns = gene.introns
        if gene.strand == "-":
            introns = mirrored(introns)
        for intron in introns:
            counter1+=1
            if gene.strand == "-":
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 11, intron[0] + 9, gene.strand, genome)
            else:
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 10, intron[0] + 10, gene.strand, genome)
            for a, base in enumerate(seq):
                if base == "A":
                    pos_matrix_5prime[0,a] = pos_matrix_5prime[0,a]+1
                if base == "C":
                    pos_matrix_5prime[1,a] = pos_matrix_5prime[1,a]+1
                if base == "T":
                    pos_matrix_5prime[2,a] = pos_matrix_5prime[2,a]+1
                if base == "G":
                    pos_matrix_5prime[3,a] = pos_matrix_5prime[3,a]+1

            if gene.strand == "-":
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[1] - 10, intron[1] + 10, gene.strand, genome)
            else:
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[1] - 11, intron[1] + 9, gene.strand, genome)
            for b, base in enumerate(seq):
                if base == "A":
                    pos_matrix_3prime[0,b] = pos_matrix_3prime[0,b]+1
                if base == "C":
                    pos_matrix_3prime[1,b] = pos_matrix_3prime[1,b]+1
                if base == "T":
                    pos_matrix_3prime[2,b] = pos_matrix_3prime[2,b]+1
                if base == "G":
                    pos_matrix_3prime[3,b] = pos_matrix_3prime[3,b]+1
    print counter1
    print counter2

    a = 0
    while a < 4:
        b = 0
        while b < 20:
            pos_matrix_5prime[a,b] = (pos_matrix_5prime[a,b])/36855.
            pos_matrix_3prime[a,b] = (pos_matrix_3prime[a,b])/36855.
            b += 1
        a += 1

    print sum(pos_matrix_5prime)        
    print pos_matrix_5prime
    print sum(pos_matrix_3prime)
    print pos_matrix_3prime
    
    return (pos_matrix_5prime, pos_matrix_3prime)


def score_new_sites(df, pos_matrix_5prime, pos_matrix_3prime):
    #Compare each splice site and add scores to dataframe - needs columns named 'coord_1', 'coord_2', 'sequence1', 'sequence2'
    df['5p score'] = ''
    df['3p score'] = ''
    df['intron length'] = 'NA'

    n = 0
    for n in range(len(df)):
        strand = df['strand'][n]
        if strand == '+':
            coord1 = int(df['coord_1'][n])
            coord2 = int(df['coord_2'][n])
            seq5 = df['sequence1'][n]
            seq3 = df['sequence2'][n]
        elif strand == '-':
            coord1 = int(df['coord_2'][n])
            coord2 = int(df['coord_1'][n])
            seq5 = df['sequence1'][n]
            seq3 = df['sequence2'][n]
        
        gene_matrix_5prime = np.zeros([4,20])
        for a, base in enumerate(seq5):
            if base == "A":
                gene_matrix_5prime[0,a] = gene_matrix_5prime[0,a]+1
            if base == "C":
                gene_matrix_5prime[1,a] = gene_matrix_5prime[1,a]+1
            if base == "T":
                gene_matrix_5prime[2,a] = gene_matrix_5prime[2,a]+1
            if base == "G":
                gene_matrix_5prime[3,a] = gene_matrix_5prime[3,a]+1

        gene_matrix_3prime = np.zeros([4,20])
        for b, base in enumerate(seq3):
            if base == "A":
                gene_matrix_3prime[0,b] = gene_matrix_3prime[0,b]+1
            if base == "C":
                gene_matrix_3prime[1,b] = gene_matrix_3prime[1,b]+1
            if base == "T":
                gene_matrix_3prime[2,b] = gene_matrix_3prime[2,b]+1
            if base == "G":
                gene_matrix_3prime[3,b] = gene_matrix_3prime[3,b]+1

        #Calculate Scores (score of 0 is perfect, higher score is worse, 40 is highest possible)
        score_5prime = 0
        score_3prime = 0
        a = 0
        while a < 4:
            b = 0
            while b < 20:
                score_5prime += abs(pos_matrix_5prime[a,b] - (gene_matrix_5prime[a,b]))
                score_3prime += abs(pos_matrix_3prime[a,b] - (gene_matrix_3prime[a,b]))
                b += 1
            a += 1

        intron_length = abs(coord2-coord1)
    
        df.loc[n,'5p score'] = score_5prime
        df.loc[n,'3p score'] = score_3prime
        df.loc[n,'intron length'] = intron_length
        
    return df

def reformat_df(df, sample_list):
    new_columns1=['Gene', 'as_event_type', 'chr', 'coord_1', 'coord_2', 'sequence1', 'sequence2', 'intron length', 'intron sequence', '5p score', '3p score', 'strand', 'intron', '#Contains_Novel_or_Only_Known(Annotated)_Junctions', 'contained in']
    new_df1 = df[new_columns1]
    new_df1 = new_df1.sort_values('as_event_type', axis=0)
    new_df1 = new_df1.rename(columns = {'as_event_type':'event','sequence1':'5p sequence','sequence2':'3p sequence','#Contains_Novel_or_Only_Known(Annotated)_Junctions':'known or novel'})
    new_df1 = new_df1.reset_index()
    new_df1 = new_df1.drop('index', axis=1)
    new_df1[['5p score','3p score']] = new_df1[['5p score','3p score']].astype(float)
    
    new_columns2 = ['Gene', 'as_event_type', 'chr', 'coord_1', 'coord_2', 'strand']+sample_list
    new_df2 = df[new_columns2]
    new_df2 = new_df2.sort_values('as_event_type', axis=0)
    new_df2 = new_df2.rename(columns = {'as_event_type':'event'})
    new_df2 = new_df2.reset_index()
    new_df2 = new_df2.drop('index', axis=1)
    new_df2[sample_list] = new_df2[sample_list].astype(float)
    
    return new_df1, new_df2