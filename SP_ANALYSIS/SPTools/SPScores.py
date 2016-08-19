'''Script for calling 5' and 3' splice sites and relative strengths of each sites 

Usage: python Call_splice_site_motifs.py

This will generate an output spreadsheet with relative splice-site strength for each transcript in the genome'''

import sys
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/')
import GeneUtility
sys.path.insert(0, '/Users/jordanburke/RNA-is-awesome/SP_ANALYSIS/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/SP_ANALYSIS/')
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
    
    junc_df['contained in'] = None
    n_max = len(junc_df.columns)
    n = 11
    while n < n_max:
        if junc_df.columns[n].startswith('avg'): continue
        elif junc_df.columns[n] == 'contained in': n += 1
        else:
            junc_df['avg_'+junc_df.columns[n]+'+'+junc_df.columns[n+1]] = junc_df[[(junc_df.columns[n]),(junc_df.columns[n+1])]].mean(axis=1)
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
        if flag == True:
            filt_df.iloc[[n]] = junc_df.iloc[[n]]
            filt_df.loc[n,'contained in'] = sample_str

    filt_df.replace(to_replace=np.NaN, value='NA', inplace=True)
    filt_df = filt_df[filt_df['strand'] != 'NA']
    
    filt_df['coord_1'] = filt_df['exclusion_junctions'].str.split(':').str.get(1).str.split('-').str.get(0)
    filt_df['coord_2'] = filt_df['exclusion_junctions'].str.split(':').str.get(1).str.split('-').str.get(1).str.split(';').str.get(0)
    
    filt_df = filt_df.reset_index()
    return filt_df

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
    df['intron'] = "Novel"
    df['sequence1'] = ''
    df['sequence2'] = ''

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
                    if coord1-2 in splice_dict[transcript][0]:
                        m = splice_dict[transcript][0].index(coord1-2)
                        if coord2-1 in splice_dict[transcript][1]:
                            if m+1 == 1:
                                df.loc[n,'intron'] = 'First'
                            elif m+1 == len(splice_dict[transcript][0]):
                                df.loc[n,'intron'] = 'Last'
                            elif m+1 > 1 and m+1 < len(splice_dict[transcript][0]):
                                df.loc[n,'intron'] = 'Middle'
                
                elif strand == '-':
                    if coord1 in splice_dict[transcript][1]:
                        m = splice_dict[transcript][1].index(coord1)
                        if coord2 in splice_dict[transcript][0]:
                            if len(splice_dict[transcript][0])-m == 1:
                                df.loc[n,'intron'] = 'First'
                            elif m == 0:
                                df.loc[n,'intron'] = 'Last'
                            elif m > 0 and m < len(splice_dict[transcript][0]):
                                df.loc[n,'intron'] = 'Middle'

        if strand == '+':
            sequence1 = fasta_dict[chrom][(coord1-3):(coord1+5)]
            sequence2 = fasta_dict[chrom][(coord2-6):(coord2+2)]
        elif strand == '-':
            sequence1 = fasta_dict[chrom][(coord2-6):(coord2+2)]
            sequence1 = SP.reverse_complement(sequence1)
            sequence2 = fasta_dict[chrom][(coord1-3):(coord1+5)]
            sequence2 = SP.reverse_complement(sequence2)
        
        df.loc[n,'sequence1'] = sequence1
        df.loc[n,'sequence2'] = sequence2

    df = df.reset_index()
    return df

def mirrored(listtree):
    return [(mirrored(x) if isinstance(x, list) else x) for x in reversed(listtree)]

def generate_consensus_matrix():
    #Populate gene dictionary and build genome
    gene_list_as_dict = GeneUtility.GeneDict()
    gene_list_as_dict.PopulateFromFile_new()
    genome = GeneUtility.BuildGenome()
    print genome.keys()

    #First generate a consensus matrix for 5' and 3' splice site, where 1st row is A counts, second row is C, third row is T, fourth row is G.
    pos_matrix_5prime = np.zeros([4,8])
    pos_matrix_3prime = np.zeros([4,8])

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
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 7, intron[0] + 1, gene.strand, genome)
            else:
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[0] - 2, intron[0] + 6, gene.strand, genome)
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
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[1] - 2, intron[1] + 6, gene.strand, genome)
            else:
                seq = GeneUtility.SequenceByGenLoc("chr{0}".format(gene.chromosome), intron[1] - 7, intron[1] + 1, gene.strand, genome)
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
        while b < 8:
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
        
        gene_matrix_5prime = np.zeros([4,8])
        
        gene_matrix_5prime = np.zeros([4,8])
        for a, base in enumerate(seq5):
            if base == "A":
                gene_matrix_5prime[0,a] = gene_matrix_5prime[0,a]+1
            if base == "C":
                gene_matrix_5prime[1,a] = gene_matrix_5prime[1,a]+1
            if base == "T":
                gene_matrix_5prime[2,a] = gene_matrix_5prime[2,a]+1
            if base == "G":
                gene_matrix_5prime[3,a] = gene_matrix_5prime[3,a]+1

        gene_matrix_3prime = np.zeros([4,8])
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
            while b < 8:
                score_5prime += abs(pos_matrix_5prime[a,b] - (gene_matrix_5prime[a,b]))
                score_3prime += abs(pos_matrix_3prime[a,b] - (gene_matrix_3prime[a,b]))
                b += 1
            a += 1

        intron_length = abs(coord2-coord1)
    
        df.loc[n,'5p score'] = score_5prime
        df.loc[n,'3p score'] = score_3prime
        df.loc[n,'intron length'] = intron_length
        
    return df

def reformat_df(df):
    new_columns=['Gene', 'as_event_type', 'chr', 'coord_1', 'coord_2', 'sequence1', 'sequence2', 'intron length', '5p score', '3p score', 'strand', 'intron', '#Contains_Novel_or_Only_Known(Annotated)_Junctions', 'contained in']
    new_df = df[new_columns]
    new_df = new_df.sort_values('as_event_type', axis=0)
    new_df = new_df.rename(columns = {'as_event_type':'event','sequence1':'5p sequence','sequence2':'3p sequence','#Contains_Novel_or_Only_Known(Annotated)_Junctions':'known or novel'})
    new_df = new_df.reset_index()
    new_df = new_df.drop('index', axis=1)
    new_df[['5p score','3p score']] = new_df[['5p score','3p score']].astype(float)
    return new_df