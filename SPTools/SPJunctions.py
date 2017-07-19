import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome')
sys.path.append('/home/jordan/RNA-is-awesome')
import SPTools as SP
import pandas as pd
import numpy as np
import pysam

def convert_chrom(chrom):
    rom_lat = {'I':'chr1','II':'chr2','III':'chr3'}
    if chrom in rom_lat:
        chrom = rom_lat[chrom]
    return chrom

def build_junction_dict(junction_bed, gff3_file, transcript_dict, organism=None):
    junction_dict = {}
    transcript_by_chr = {}
    unassigned_count = 0
    
    ss_dict, flag = SP.list_splice_sites(gff3_file, organism=organism)   
    ss_by_gene = SP.collapse_ss_dict(ss_dict)
    ss_by_gene = {k:v for k, v in ss_by_gene.items() if len(v) > 0}
    
    for transcript, coords in transcript_dict.iteritems():
        if transcript[:-2] in ss_by_gene.keys():
            chromosome = coords[3]
            #junction_dict[transcript] = []
            if chromosome in transcript_by_chr:
                transcript_by_chr[chromosome].append(transcript)
            else:
                transcript_by_chr[chromosome] = []
                transcript_by_chr[chromosome].append(transcript)

    a = 0
    with open(junction_bed, 'r') as fin:
        for line in fin:
            a += 1
            jct_transcript = None
            jct_type = 'Other'
            intron_num = None
            if line.startswith('c') or line.startswith('I'):
                columns = line.split('\t')
                lat_rom = {'I':'chr1', 'II':'chr2', 'III':'chr3'}
                chromosome = columns[0]
                if chromosome in lat_rom:
                    chromosome = lat_rom[chromosome]
                
                if chromosome in transcript_by_chr:
                    transcript_list = transcript_by_chr[chromosome]
                    #if organism == 'pombe':
                    #    transcript_list = [x for x in transcript_list if x[:-2] in ss_by_gene.keys()]
                
                strand = columns[5]
                if strand == '+':
                    jct_start = int(columns[1])+int(columns[10].split(',')[0])-1
                    jct_end = int(columns[2])-int(columns[10].split(',')[1])-1
                elif strand == '-':
                    jct_start = int(columns[2])-int(columns[10].split(',')[1])
                    jct_end = int(columns[1])+int(columns[10].split(',')[0])
                depth = int(columns[4])
                size = abs(jct_end-jct_start)
                
                assigned = False
                for transcript in transcript_list:
                    if jct_start > transcript_dict[transcript][0] and jct_end < transcript_dict[transcript][1] and strand == transcript_dict[transcript][2]:
                        assigned = True
                        jct_transcript = transcript
                        all_sites = zip(*ss_by_gene[transcript[:-2]])
                        try:
                            if jct_start in all_sites[0] and jct_end in all_sites[1]:
                                jct_type = 'Annotated'
                                ann_size = size
                                intron_num = all_sites[0].index(jct_start)+1
                                ann_start = jct_start
                                ann_stop = jct_end
                                break
                            else:
                                n=0
                                for intron in ss_by_gene[transcript[:-2]]:
                                    n += 1
                                    ann_size = None
                                    if strand == '+':
                                        if jct_start > intron[0] and jct_end < intron[1]:
                                            ann_size = abs(intron[1]-intron[0])
                                            jct_type = 'Nested'
                                            ann_start = intron[0]
                                            ann_stop = intron[1]
                                            intron_num = n
                                            break
                                        elif jct_start >= intron[0] and jct_end == intron[1]:
                                            ann_size = abs(intron[1]-intron[0])
                                            jct_type = '3p tethered'
                                            ann_start = intron[0]
                                            ann_stop = intron[1]
                                            intron_num = n
                                            break
                                        elif jct_start == intron[0] and jct_end <= intron[1]:
                                            ann_size = abs(intron[1]-intron[0])
                                            jct_type = '5p tethered'
                                            ann_start = intron[0]
                                            ann_stop = intron[1]
                                            intron_num = n
                                            break
                                        
                                    elif strand == '-':
                                        if jct_start < intron[0] and jct_end > intron[1]:
                                            jct_type = 'Nested'
                                            ann_size = intron[0]-intron[1]
                                            ann_start = intron[0]
                                            ann_stop = intron[1]
                                            intron_num = n
                                            break
                                        elif jct_start <= intron[0] and jct_end == intron[1]:
                                            jct_type = '5p tethered'
                                            ann_size = intron[0]-intron[1]
                                            ann_start = intron[0]
                                            ann_stop = intron[1]
                                            intron_num = n
                                            break
                                        elif jct_start == intron[0] and jct_end >= intron[1]:
                                            jct_type = '3p tethered'
                                            ann_size = intron[0]-intron[1]
                                            ann_start = intron[0]
                                            ann_stop = intron[1]
                                            intron_num = n
                                            break
                            break
                        except IndexError:
                            print transcript
                if assigned is False: unassigned_count += 1
 
                try:
                    if jct_transcript != None:
                        if ann_size == None:
                            jct_type = "Other"
                            ann_size = 0
                            ann_start = None
                            ann_stop = None
                            intron_num = None
                        if (jct_transcript, ann_size) not in junction_dict:
                            junction_dict[(jct_transcript, ann_size)] = []
                        junction_dict[(jct_transcript, ann_size)].append([chromosome, jct_start, jct_end, strand, depth, jct_type, size, ann_size, ann_start, ann_stop])
                except ValueError:
                    print jct_transcript
                    print jct_type

    print str(unassigned_count)+' junctions not assigned to transcripts'
    return junction_dict

def build_junction_df(junction_bed, gff3_file, fasta, organism=None):
    transcript_dict = SP.build_transcript_dict(gff3_file, organism=organism)
    if type(fasta) == str:
        fasta=SP.make_fasta_dict(fasta)
    junction_dict = build_junction_dict(junction_bed, gff3_file, transcript_dict, organism=organism)
    junction_count = 0
    for tx, junctions in junction_dict.iteritems():
        junction_count += len(junctions)
    
    junction_df = pd.DataFrame(index=range(junction_count), columns=['intron tuple','chromosome','start','end','strand','depth','type','size','annotated intron size','annotated intron start','annotated intron end'])
    n=0
    for tx, junctions in junction_dict.iteritems():
        for junction in junctions:
            junction_df.ix[n] = [tx]+junction
            n+=1
    
    sequence1 = []
    sequence2 = []
    ann_seq1 = []
    ann_seq2 = []
    seq_type1 = []
    seq_type2 = []
    df_tx = []
    for index, row in junction_df.iterrows():
        df_tx.append(row['intron tuple'][0])
        chrom = convert_chrom(row['chromosome'])
        if row['strand'] == '+':
            curr1 = fasta[chrom][(row['start']-1):(row['start']+7)]
            sequence1.append(curr1)
            curr2 = fasta[chrom][(row['end']-5):(row['end']+3)]
            sequence2.append(curr2)
            if row['annotated intron start'] is None:
                ann_seq1.append(None)
                ann_seq2.append(None)
            else:
                ann_seq1.append(fasta[chrom][(row['annotated intron start']-1):(row['annotated intron start']+7)])
                ann_seq2.append(fasta[chrom][(row['annotated intron end']-5):(row['annotated intron end']+3)])
        elif row['strand'] == '-':
            curr1 = SP.reverse_complement(fasta[chrom][(row['start']-6):(row['start']+2)])
            sequence1.append(curr1)
            curr2 = SP.reverse_complement(fasta[chrom][(row['end']-2):(row['end']+6)])
            sequence2.append(curr2)
            if row['annotated intron start'] is None:
                ann_seq1.append(None)
                ann_seq2.append(None)
            else:
                ann_seq1.append(SP.reverse_complement(fasta[chrom][row['annotated intron start']-6:row['annotated intron start']+2]))
                ann_seq2.append(SP.reverse_complement(fasta[chrom][row['annotated intron end']-2:row['annotated intron end']+6]))
        else:
            sequence1.append('NNNNNNNN')
            sequence2.append('NNNNNNNN')
            ann_seq1.append('NNNNNNNN')
            ann_seq2.append('NNNNNNNN')
        
        
        if row['type'] == 'Annotated': 
            seq_type1.append('5p annotated')
            seq_type2.append('3p annotated')
        elif row['type'] == '5p tethered':
            seq_type1.append('5p annotated')
            seq_type2.append(curr2[4:6])
        else:
            seq_type1.append(curr1[2:4])
            seq_type2.append(curr2[4:6])
            
    junc_seq_df = junction_df
    junc_seq_df['sequence1'] = sequence1
    junc_seq_df['sequence2'] = sequence2
    junc_seq_df['seq type1'] = seq_type1
    junc_seq_df['seq type2'] = seq_type2
    junc_seq_df['annotated sequence1'] = ann_seq1
    junc_seq_df['annotated sequence2'] = ann_seq2
    junc_seq_df['transcript'] = df_tx
    
    return junc_seq_df
            
def combine_junctions(junc_df1, junc_df2):
    junc_df1['genome coord'] = junc_df1['chromosome'].str.cat(junc_df1['start'].values.astype(str), sep=':').str.cat(junc_df1['end'].values.astype(str), sep='-')
    junc_df2['genome coord'] = junc_df2['chromosome'].str.cat(junc_df2['start'].values.astype(str), sep=':').str.cat(junc_df2['end'].values.astype(str), sep='-')
    new_df = junc_df1.append(junc_df2[~junc_df2['genome coord'].isin(junc_df1['genome coord'].tolist())])
    new_df.drop('genome coord', inplace=True)
    return new_df
    

def add_int_levels(junction_df, int_levels_df, gff3, organism=None):
    ss_dict, flag = SP.list_splice_sites(gff3, organism=organism)
    tx_dict = SP.build_transcript_dict(gff3, organism=organism)
    intron_sizes = {}
    for transcript, sites in ss_dict.iteritems():
        if len(sites[0]) > 0:
            intron_sizes[transcript] = []
            if tx_dict[transcript][2] == '+':
                n=0
                for n in range(len(sites[0])):
                    intron_sizes[transcript].append((n+1,sites[1][n]-sites[0][n]))
            elif tx_dict[transcript][2] == '-':
                n=0
                for n in range(len(sites[0])):
                    intron_sizes[transcript].append((len(sites[0])-n,sites[0][n]-sites[1][n]))
    
    int_level_columns = {}
    for entry in int_levels_df.columns:
        if entry[1] == '5prime Normalized':
            int_level_columns[entry] = []
    
    for index, row in junction_df.iterrows():
        tx = row['intron tuple'][0]
        intron_size = row['intron tuple'][1]
        if intron_size == 0:
            for sample in int_level_columns:
                int_level_columns[sample].append(np.NaN)
        else:
            try:
                intron_ix = zip(*intron_sizes[tx])[1].index(intron_size)
                intron_num = zip(*intron_sizes[tx])[0][intron_ix]
                
                for sample in int_level_columns:
                    int_level_columns[sample].append(int_levels_df.loc[(tx,intron_num)][sample])
                    
            except ValueError:
                for sample in int_level_columns:
                    int_level_columns[sample].append(np.NaN)
            except KeyError:
                for sample in int_level_columns:
                    int_level_columns[sample].append(np.NaN)
    for sample in int_level_columns:
        junction_df[sample] = int_level_columns[sample]
    
    junction_df = junction_df.dropna(how='any')
    return junction_df

def add_occupancy(junction_df, occupancy_df):
    occ_columns = {}
    for entry in occupancy_df.columns:
        if entry[1] == 'Normalized to mature':
            occ_columns[entry] = []
    
    for index, row in junction_df.iterrows():
        tx = row['transcript']
        try:
            for sample in occ_columns:
                occ_columns[sample].append(occupancy_df.loc[tx][sample])

        #except ValueError:
        #    for sample in int_level_columns:
        #        int_level_columns[sample].append(np.NaN)
        except KeyError:
            for sample in occ_columns:
                occ_columns[sample].append(np.NaN)
    
    for sample in occ_columns:
        junction_df[sample] = occ_columns[sample]
        
    junction_df = junction_df.dropna(how='any')
    return junction_df