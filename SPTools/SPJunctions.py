import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/SPTools')
import SPTables
import SPPeaks
import SPPlots
import SPScores
import pandas as pd
import numpy as np
import pysam

def convert_chrom(chrom):
    rom_lat = {'I':'chr1','II':'chr2','III':'chr3'}
    if chrom in rom_lat:
        chrom = rom_lat[chrom]
    return chrom

def collapse_ss_dict(splice_site_dict):
    ss_by_gene = {}
    transcripts = splice_site_dict.keys()
    transcripts = [x[:-2] for x in transcripts]
    genes = list(set(transcripts))
    for gene in genes:
        for transcript, sites in splice_site_dict.iteritems():
            if gene in transcript:
                if gene not in ss_by_gene:
                    ss_by_gene[gene] = set()
                ss_by_gene[gene].update(set(zip(sites[0],sites[1])))
    return ss_by_gene

def build_junction_dict(junction_bed, gff3_file, transcript_dict, organism=None):
    junction_dict = {}
    transcript_by_chr = {}
    unassigned_count = 0
    
    ss_dict, flag = SPPeaks.list_splice_sites(gff3_file, organism=organism)   
    ss_by_gene = collapse_ss_dict(ss_dict)
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
            if a%1000 == 0:
                print a
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
    transcript_dict = SPPeaks.build_transcript_dict(gff3_file, organism=organism)
    if type(fasta) == str:
        fasta=SPScores.make_fasta_dict(fasta)
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
            curr1 = SPPeaks.reverse_complement(fasta[chrom][(row['start']-6):(row['start']+2)])
            sequence1.append(curr1)
            curr2 = SPPeaks.reverse_complement(fasta[chrom][(row['end']-2):(row['end']+6)])
            sequence2.append(curr2)
            if row['annotated intron start'] is None:
                ann_seq1.append(None)
                ann_seq2.append(None)
            else:
                ann_seq1.append(SPPeaks.reverse_complement(fasta[chrom][row['annotated intron start']-6:row['annotated intron start']+2]))
                ann_seq2.append(SPPeaks.reverse_complement(fasta[chrom][row['annotated intron end']-2:row['annotated intron end']+6]))
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
    print len(junc_df1)
    print len(junc_df2)
    junc_df1['genome coord'] = junc_df1['chromosome'].str.cat(junc_df1['start'].values.astype(str), sep=':').str.cat(junc_df1['end'].values.astype(str), sep='-')
    junc_df2['genome coord'] = junc_df2['chromosome'].str.cat(junc_df2['start'].values.astype(str), sep=':').str.cat(junc_df2['end'].values.astype(str), sep='-')
    new_df = junc_df1.append(junc_df2[~junc_df2['genome coord'].isin(junc_df1['genome coord'].tolist())])
    new_df.drop('genome coord', inplace=True)
    print len(new_df)
    return new_df
    
    
def junction_v_peak(junction_dict, peak_pipeline_out):
    new_dict = {}
    peaks_df = pd.read_csv(peak_pipeline_out, sep='\t')
    peaks_df.drop_duplicates(inplace=True)
    #print peaks_df
    novel_peak_dict = {}
    junction_dict = {k:v for k, v in junction_dict.items() if type(k) == tuple}
    for intron in junction_dict.keys():
        transcript_df = peaks_df[peaks_df['transcript'] == intron[0][:-2]]
        peaks = transcript_df[' coordinate'].tolist()
        peaks = map(int, peaks)
        n = 0
        matched_peaks = set()
        for junction in junction_dict[intron]:
            peak_at_jct5 = False
            peak_at_jct3 = False

            if junction[0]+1 in peaks:
                peak_at_jct5 = True
                matched_peaks.add(junction[0])
            if junction[1]+1 in peaks:
                peak_at_jct3 = True
                matched_peaks.add(junction[1])
            
            if intron not in new_dict:
                new_dict[intron] = []
            new_dict[intron].append(junction_dict[intron][n]+[peak_at_jct5,peak_at_jct3])
            n += 1
        unmatched_peaks = list(set(peaks).difference(matched_peaks))
        
        if intron[0] not in novel_peak_dict:
            novel_peak_dict[intron[0]] = []
        for row in transcript_df.iterrows():
            if row[1][' coordinate'] in unmatched_peaks:
                unmatched_peaks.remove(row[1][' coordinate'])
                novel_peak_dict[intron[0]].append(list(row[1]))
    return new_dict, novel_peak_dict

def count_junctions(jct_peak_dict, novel_peaks):
    jct_counter = 0
    ann_counter = 0
    unk_counter = 0
    peak_at_ann = 0
    peak_at_unk = 0
    rec_count = 0
    transcript_count = 0
    total_rec = 0
    
    rec_dict = {}
    ann_dict = {}
    for intron, jcts in jct_peak_dict.iteritems():
        recursive = False
        rec_per_t = 0
        if len(jcts) > 0: transcript_count += 1
        sizes = []
        reads = []
        jct_sizes = []
        ann_sizes = []
        intron_num = []
        for jct in jcts:
            jct_counter += 1
            if jct[3] == 'Annotated': 
                ann_counter += 1
                if jct[2] > 10:
                    sizes.append(jct[4])
                    reads.append(jct[2])
                if jct[9] == True or jct[10] == True: peak_at_ann += 1
            else:
                unk_counter += 1
                if jct[9] == True or jct[10] == True: peak_at_unk += 1
                    
            if jct[3] == 'Nested' or jct[3] == '5p tethered' or jct[3] == '3p tethered':
                recursive = True
                rec_per_t += 1
                jct_sizes.append(jct[4])
                ann_sizes.append(jct[5])
                intron_num.append(jct[6])
        ann_dict[intron] = [sizes, reads]
        if recursive == True: 
            rec_count += 1
            rec_dict[intron] = [rec_per_t,jct_sizes,ann_sizes, intron_num]
            total_rec += rec_per_t
    print str(jct_counter)+' junctions detected'
    print str(ann_counter)+' junctions at annotated sites'
    print str(unk_counter)+' junctions at novel sites'    
    print str(peak_at_ann)+" junctions corresponding to peaks at annotated sites"
    print str(peak_at_unk)+" junctions corresponding to peaks at novel sites\n"
    
    print str(rec_count)+" introns with recursive events out of "+str(transcript_count)
    print str(total_rec)+" total recursive events detected\n"

    
    peak_counter = 0
    ann_peak_counter = 0
    novel_peak_counter = 0
    GU = 0
    AG = 0
    for transcript, peaks in novel_peaks.iteritems():
        for peak in peaks:
            peak_counter += 1
            if peak[2] != 'Unknown': ann_peak_counter += 1
            else:
                novel_peak_counter += 1
                if peak[8] == "5'": GU += 1
                elif peak[8] == "3'": AG += 1
    #print str(peak_counter)+' peaks not at junctions'
    #print 'Of these:'
    #print str(ann_peak_counter)+' peaks at annotated sites'
    #print str(novel_peak_counter)+' unpredicted peaks'
    #print str(GU)+' unpredicted peaks with 5prime splice sites'
    #print str(AG)+' unpredicted peaks with 3prime splice sites'
    
    return rec_dict, ann_dict

def rec_seq(jct_dict, gff3_file, fasta_file):
    rec_dict = {}
    for intron, junctions in jct_dict.iteritems():
        transcript = intron[0]
        for junction in junctions:
            if transcript not in rec_dict:
                rec_dict[transcript] = [[],[]]
            rec_dict[transcript][0].append(junction[0])
            rec_dict[transcript][1].append(junction[1])
    seq_dict = SPScores.get_sequence(rec_dict, gff3_file, fasta_file)
    return seq_dict

def read_kallisto_abundance(abundance_file):
    kallisto_dict = {}
    #kallisto_dict[transcript] = [length, eff_length, est_counts, tpm]
    with open(abundance_file, 'r') as fin:
        for line in fin:
            if line.startswith('target'):
                continue
            columns = line.split('\t')
            kallisto_dict[columns[0]] = [int(columns[1]), float(columns[2]), float(columns[3]), float(columns[4].strip())]
    return kallisto_dict

def compare_reads_at_junctions(bam_file, gff3, transcript_list='All'):
    transcript_dict = SPPeaks.build_transcript_dict(gff3)
    ss_dict, flag = SPPeaks.list_splice_sites(gff3)
    if transcript_list != 'All':
        transcript_dict = {k:v for k,v in transcript_dict.items() if k in transcript_list}
    
    reads_df = pd.DataFrame(columns=['5p splice site reads', '3p splice site reads', 'Exon junction reads', 'Reads in transcript', 'Transcript length'], index=transcript_dict.keys())
    bam = pysam.Samfile(bam_file)
    for tx, info in transcript_dict.iteritems():
        print tx
        chrom = transcript_dict[tx][3]
        start = transcript_dict[tx][0]
        end = transcript_dict[tx][1]
        strand = transcript_dict[tx][2]
        length = abs(end-start)
        
        tx_reads = 0
        five_reads = 0
        three_reads = 0 
        e_e_reads = 0
        tx_iter = bam.fetch(chrom, start, end)
        for read in tx_iter:
            if read.is_reverse and strand == '+':
                tx_reads += 1
            elif not read.is_reverse and strand == '-':
                tx_reads += 1                    
                            
            if 'N' not in read.cigarstring:
                for site in ss_dict[tx][0]:
                    if site in range(read.reference_start,read.reference_end):
                        if read.is_reverse and strand == '+':
                            five_reads += 1
                        elif not read.is_reverse and strand == '-':
                            five_reads += 1
                for site in ss_dict[tx][1]:
                    if site in range(read.reference_start, read.reference_end):
                        if read.is_reverse and strand == '+':
                            three_reads += 1
                        elif not read.is_reverse and strand == '-':
                            three_reads += 1
            elif 'N' in read.cigarstring:
                e_e_reads += 1
            
            row = [five_reads, three_reads, e_e_reads, tx_reads, length]
            reads_df.loc[tx] = row
    return reads_df

def count_sense_vs_antisense(bam_file, gff3, transcript_list='All'):
    transcript_dict = SPPeaks.build_transcript_dict(gff3)
    if transcript_list != 'All':
        transcript_dict = {k:v for k,v in transcript_dict.items() if k in transcript_list}
        
    reads_df = pd.DataFrame(columns=['Sense reads', 'Antisense reads', 'Transcript length'], index=transcript_dict.keys())
    bam = pysam.Samfile(bam_file)
    for tx, info in transcript_dict.iteritems():
        print tx
        chrom = transcript_dict[tx][3]
        start = transcript_dict[tx][0]
        end = transcript_dict[tx][1]
        strand = transcript_dict[tx][2]
        length = abs(end-start)
        
        sense_reads = 0
        anti_reads = 0
        tx_iter = bam.fetch(chrom, start, end)
        for read in tx_iter:
            bool_tuple = (read.is_reverse, read.is_read1)
            if read.is_reverse and read.is_read1:
                if strand == '+':
                    sense_reads += 1
                elif strand == '-':
                    anti_reads += 1
            elif read.is_reverse and read.is_read2:
                if strand == '-':
                    sense_reads += 1
                elif strand == '+':
                    anti_reads += 1
            elif not read.is_reverse and read.is_read1:
                if strand == '-':
                    sense_reads += 1
                elif strand == '+':
                    anti_reads += 1
            elif not read.is_reverse and read.is_read2:
                if strand == '+':
                    sense_reads += 1
                elif strand == '-':
                    anti_reads += 1
                    
            row = [sense_reads, anti_reads, length]
            reads_df.loc[tx] = row
    return reads_df

def make_seq_df(jct_dict, transcript_dict, fasta_dict):
    jct_count = 0
    for intron, junctions in jct_dict.iteritems():
        if intron[1] > 0:
            for junction in junctions:
                if junction[3] != 'other':
                    jct_count += 1
    columns = ['transcript','annotated intron size','chromosome','coord_1','coord_2','strand','size','junction type','sequence1','sequence2']
    df = pd.DataFrame(columns=columns, index= range(jct_count))
    df_ann = pd.DataFrame(columns=columns, index=range(jct_count))
    print df.shape
    n = 0
    for intron, junctions in jct_dict.iteritems():
        if type(intron) != tuple: 
            pass
        else:
            transcript = intron[0]
            ann_size = intron[1]
            if ann_size > 0:
                chrom = transcript_dict[transcript][3]
                strand = transcript_dict[transcript][2]

                for junction in junctions:
                    ann_start = junction[6]
                    ann_stop = junction[7]
                    coord_1 = junction[0]
                    coord_2 = junction[1]
                    size = junction[4]
                    jct_type = junction[3]

                    if jct_type != 'other':  
                        if strand == '+':
                            sequence1 = fasta_dict[chrom][(coord_1-1):(coord_1+7)]
                            sequence2 = fasta_dict[chrom][(coord_2-5):(coord_2+3)]
                            ann_seq1 = fasta_dict[chrom][(ann_start-1):(ann_start+7)]
                            ann_seq2 = fasta_dict[chrom][(ann_stop-5):(ann_stop+3)]

                        elif strand == '-':
                            sequence1 = fasta_dict[chrom][(coord_1-6):(coord_1+2)]
                            sequence1 = SPPeaks.reverse_complement(sequence1)
                            sequence2 = fasta_dict[chrom][(coord_2-2):(coord_2+6)]
                            sequence2 = SPPeaks.reverse_complement(sequence2)

                            ann_seq1 = fasta_dict[chrom][ann_start-6:ann_start+2]
                            ann_seq2 = fasta_dict[chrom][ann_stop-2:ann_stop+6]
                            ann_seq1 = SPPeaks.reverse_complement(ann_seq1)
                            ann_seq2 = SPPeaks.reverse_complement(ann_seq2)

                        row = [transcript, ann_size, chrom, coord_1, coord_2, strand, size, jct_type, sequence1, sequence2]
                        df.iloc[n] = row

                        ann_row = [transcript, ann_size, chrom, ann_start, ann_stop, strand, ann_size, jct_type, ann_seq1, ann_seq2]
                        df_ann.iloc[n] = ann_row

                        n += 1
    
    df_ann = df_ann.drop_duplicates()
    df_ann = df_ann.reset_index()
    df_ann = df_ann.drop('index', axis=1)
    return df, df_ann

def add_int_levels(junction_df, int_levels_df, gff3, organism=None):
    ss_dict, flag = SPPeaks.list_splice_sites(gff3, organism=organism)
    tx_dict = SPPeaks.build_transcript_dict(gff3, organism=organism)
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
            #print entry
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