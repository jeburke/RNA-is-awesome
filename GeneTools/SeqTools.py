import sys
import os
from subprocess import check_output
import math
import numpy as np
import pandas as pd
sys.path.insert(0, '/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.insert(0, '/home/jordan/RNA-is-awesome/')
sys.path.insert(0, '/Users/jordanburke/CodeBase/RNA-is-awesome/')
import GeneUtility
import SPTools as SP
from collections import OrderedDict
import csv
import json

def build_transcript_dict(gff3_file, organism=None):
    with open(gff3_file,"r") as gff3:
        transcript_dict = {}
        for line in gff3:
            columns = line.split('\t')
            
            if organism == 'pombe' and len(columns) > 1:
                chr_rom = columns[0]
                rom_lat = {'I':'chr1','II':'chr2','III':'chr3','MT':'MT'}
                chrom = rom_lat[chr_rom]
                transcript_types = ['transcript','pseudogene','rRNA','snoRNA','tRNA','snRNA']
                if columns[2] in transcript_types:
                    if columns[8].split(':')[0].split('=')[1] == 'gene': continue
                    transcript = columns[8].split(';')[0].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], chrom, [], []]
                elif columns[2] == 'CDS':
                    transcript = columns[8].split(':')[1]
                    #if transcript[-2] != 'T': transcript = transcript[:-1]+'T1'
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
                        
            if len(columns) == 9 and organism is None:
                if columns[2] == "mRNA" or columns[2] == "snoRNA_gene" or columns[2] == "tRNA_gene":
                    transcript = columns[8]
                    transcript = transcript.split("=")[1]
                    transcript = transcript.split(";")[0]
                    if transcript.endswith("mRNA"): transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    #Transcript dictionary: keys are transcript, values are [start, end, strand, chromosome, CDS start, CDS end]
                    if transcript not in transcript_dict:
                        transcript_dict[transcript] = [int(columns[3]), int(columns[4]), columns[6], columns[0], [], []]
                    else:
                        transcript_dict[transcript][0] = int(columns[3])
                        transcript_dict[transcript][1] = int(columns[4])
                        transcript_dict[transcript][2] = columns[6]
                        transcript_dict[transcript][3] = columns[0]
                elif columns[2] == "CDS":
                    transcript = columns[8].split("=")[1].split(".")[0].split(';')[0]
                    if 'mRNA' in transcript: transcript = transcript.split("_")[0]
                    if transcript[-2] != 'T': transcript = transcript+'T0'
                    if transcript not in transcript_dict:
                        strand = columns[6]
                        chrom = columns[0]
                        transcript_dict[transcript] = [0,0,strand,chrom,[],[]]
                    transcript_dict[transcript][4].append(int(columns[3]))
                    transcript_dict[transcript][5].append(int(columns[4]))
    
    for tx in transcript_dict:
        if transcript_dict[tx][0] == 0:
            transcript_dict[tx][0] = transcript_dict[tx][4][0]
            transcript_dict[tx][1] = transcript_dict[tx][5][0]
    
    transcript_dict = OrderedDict(sorted(transcript_dict.items()))
    return transcript_dict


def seq_simple(chrom, start, end, strand, fasta_dict):
    if type(fasta_dict) == str:
        with open(fasta_dict, 'r') as f:
            fasta_dict = json.load(f)
    seq = fasta_dict[chrom][start:end+1]
    if strand == '-':
        seq = SP.reverse_complement(seq)
    return seq
            
def get_peak_sequence(input_file, fasta_file, gff3_file, window=1000):
    #File where 1st column is transcript (3P prefix sometimes), 2nd column is chr, 3rd column is peak center
    tx_dict = SP.build_transcript_dict(gff3_file)
    if type(fasta_file) == dict:
        fa_dict = fasta_file
    else:
        fa_dict = SP.make_fasta_dict(fasta_file)
    seq_list = []
    with open(input_file,'r') as csv_file:
        f = csv.reader(csv_file, dialect=csv.excel)
        for row in f:
            tx = row[0]+'T0'
            if tx.startswith('3P'): tx = tx.split('3P')[1]
            chrom = 'chr'+row[1]
            try:
                center = int(row[2])
                if tx in tx_dict:
                    strand = tx_dict[tx][2]
                    start = center-window/2
                    end = center+window/2
                    seq = seq_simple(chrom, start, end, strand, fa_dict)
                    seq_list.append((tx,seq))
                else:
                    print tx+" not in GFF3 file"
            except ValueError:
                pass
    with open('{0}_peak_sequences.fa'.format(input_file.split('/')[-1].split('.')[0]),'w') as fout:
        for tx, seq in seq_list:
            fout.write('>'+tx+'\n')
            fout.write(seq+'\n')
    return seq_list


# This function reads my sorted-by-gene bedgraph format (output of SPTools.build_bedgraph_dict). Make sure the same
# transcript_dictionary is used for both function calls. In this case I have modified the dictionary to include 300 bp
# downstrame of the stop codon for each transcript

def read_bg_dict(bedgraph_dict_output, transcript_dict):
    bg_dict = {}
    count = 0
    dtype = [('coord', int), ('height', float)]
    with open(bedgraph_dict_output,'r') as f:
        n = 0
        for line in f:
            n += 1
            if len(line) > 1:
                
                #Read the transcript line
                if n%3 == 1:
                    tx = line.strip()
                    count += 1
                    #print count
                    if tx[-2] != 'T':
                        tx = tx+'T0'
                
                #Read the coordinate line
                elif n%3 == 2:
                    coords  = map(int, line.strip().split('\t'))
                
                #Read the values line
                elif n%3 == 0:
                    heights = map(float, line.strip().split('\t'))

                    all_coords = set(range(transcript_dict[tx][0],transcript_dict[tx][1]))
                    missing = all_coords.difference(coords)
                    coords = coords + list(missing)

                    #Fill in missing coordinates with zeros
                    zero_fill = [0]*len(missing)
                    heights = heights + zero_fill
                    
                    #Create a pandas series with all coordinates and sort so zeros are inserted appropriately
                    entry = pd.Series(heights, index=coords)
                    entry.sort_index(inplace=True)

                    bg_dict[tx] = entry
    return bg_dict

def seq_file_from_df(df, column_name, file_name):
    with open(file_name,'w') as fout:
        seq_list = df[column_name].tolist()
        for seq in seq_list:
            fout.write(seq+'\n')