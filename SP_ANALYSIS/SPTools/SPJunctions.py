import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/SP_ANALYSIS/SPTools')
import SPTables
import SPPeaks
import SPPlots
import SPScores
import pandas as pd
import numpy as np

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

def build_junction_dict(junction_bed, gff3_file, transcript_dict):
    junction_dict = {}
    transcript_by_chr = {}
    for transcript, coords in transcript_dict.iteritems():
        chromosome = coords[3]
        junction_dict[transcript] = []
        if chromosome in transcript_by_chr:
            transcript_by_chr[chromosome].append(transcript)
        else:
            transcript_by_chr[chromosome] = []
            transcript_by_chr[chromosome].append(transcript)
    ss_dict, flag = SPPeaks.list_splice_sites(gff3_file)   
    ss_by_gene = collapse_ss_dict(ss_dict)
    
    with open(junction_bed, 'r') as fin:
        for line in fin:
            jct_transcript = None
            jct_type = 'Other'
            if line.startswith('c'):
                columns = line.split('\t')
                chromosome = columns[0]
                if chromosome in transcript_by_chr:
                    transcript_list = transcript_by_chr[chromosome]
                strand = columns[5]
                if strand == '+':
                    jct_start = int(columns[1])+int(columns[10].split(',')[0])-1
                    jct_end = int(columns[2])-int(columns[10].split(',')[1])-1
                elif strand == '-':
                    jct_start = int(columns[2])-int(columns[10].split(',')[1])
                    jct_end = int(columns[1])+int(columns[10].split(',')[0])
                depth = int(columns[4])
                size = abs(jct_end-jct_start)
                
                for transcript in transcript_list:
                    if jct_start > transcript_dict[transcript][0] and jct_end < transcript_dict[transcript][1] and strand == transcript_dict[transcript][2]:
                        jct_transcript = transcript
                        all_sites = zip(*ss_by_gene[transcript[:-2]])
                        try:
                            if jct_start in all_sites[0] and jct_end in all_sites[1]:
                                jct_type = 'Annotated'
                                ann_size = size
                            else:
                                for intron in ss_by_gene[transcript[:-2]]:
                                    ann_size = abs(intron[1]-intron[0])
                                    if strand == '+':
                                        if jct_start >= intron[0] and jct_end <= intron[1]:
                                            jct_type = 'Intronic'
                                            if intron[1]-intron[0] > 20:
                                                jct_type = 'Recursive'
                                    elif strand == '-':
                                        if jct_start <= intron[0] and jct_end >= intron[1]:
                                            jct_type = 'Intronic'
                                            if intron[0]-intron[1] > 20:
                                                jct_type = 'Recursive'
                            break
                        except IndexError:
                            print transcript
 
                try:
                    if jct_transcript != None:
                        if jct_transcript not in junction_dict:
                            junction_dict[jct_transcript] = []
                        else:
                            junction_dict[jct_transcript].append([jct_start, jct_end, depth, jct_type, size, ann_size])
                except ValueError:
                    print jct_transcript
                    print jct_type
                    #print jct_transcript
                    #print ss_dict[jct_transcript]
                    #print junction_dict[jct_transcript]
    return junction_dict                       

def junction_v_peak(junction_dict, peak_pipeline_out):
    new_dict = {}
    peaks_df = pd.read_csv(peak_pipeline_out, sep='\t')
    novel_peak_dict = {}
    for transcript in junction_dict:
        transcript_df = peaks_df[peaks_df['transcript'] == transcript[:-2]]
        peaks = transcript_df[' coordinate'].tolist()
        peaks = map(int, peaks)
        n = 0
        matched_peaks = set()
        for junction in junction_dict[transcript]:
            peak_at_jct5 = False
            peak_at_jct3 = False
            #if transcript == 'CNAG_00483T0': 
            #    print junction[0]
            #    print junction[1]
            
            if junction[0] in peaks:
                peak_at_jct5 = True
                matched_peaks.add(junction[0])
            if junction[1] in peaks:
                peak_at_jct3 = True
                matched_peaks.add(junction[1])
            
            if transcript not in new_dict:
                new_dict[transcript] = []
            new_dict[transcript].append(junction_dict[transcript][n]+[peak_at_jct5,peak_at_jct3])
            n += 1
        unmatched_peaks = list(set(peaks).difference(matched_peaks))
        #if transcript == 'CNAG_00483T0':
        #    print junction_dict['CNAG_00483T0']
        #    print peaks
        #    print matched_peaks
        #    print unmatched_peaks
        if transcript not in novel_peak_dict:
            novel_peak_dict[transcript] = []
        for row in transcript_df.iterrows():
            if row[1][' coordinate'] in unmatched_peaks:
                novel_peak_dict[transcript].append(list(row[1]))
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
    for transcript, jcts in jct_peak_dict.iteritems():
        recursive = False
        rec_per_t = 0
        if len(jcts) > 0: transcript_count += 1
        jct_sizes = []
        ann_sizes = []
        for jct in jcts:
            jct_counter += 1
            if jct[3] == 'Annotated': 
                ann_counter += 1
                if jct[6] == True or jct[7] == True: peak_at_ann += 1
            else:
                unk_counter += 1
                if jct[6] == True or jct[7] == True: peak_at_unk += 1
                    
            if jct[3] == 'Recursive':
                recursive = True
                rec_per_t += 1
                jct_sizes.append(jct[4])
                ann_sizes.append(jct[5])
        if recursive == True: 
            rec_count += 1
            rec_dict[transcript] = [rec_per_t,jct_sizes,ann_sizes]
            total_rec += rec_per_t
    print str(jct_counter)+' junctions detected'
    print str(ann_counter)+' junctions at annotated sites'
    print str(unk_counter)+' junctions at novel sites'    
    print str(peak_at_ann)+" junctions corresponding to peaks at annotated sites"
    print str(peak_at_unk)+" junctions corresponding to peaks at novel sites\n"
    
    print str(rec_count)+" transcripts with recursive events out of "+str(transcript_count)
    print str(total_rec)+" total recursive events detected\n"

    
    peak_counter = 0
    ann_peak_counter = 0
    GU = 0
    AG = 0
    for transcript, peaks in novel_peaks.iteritems():
        for peak in peaks:
            peak_counter += 1
            if peak[2] != 'Unknown': ann_peak_counter += 1
            elif peak[8] == "5'": GU += 1
            elif peak[8] == "3'": AG += 1
    print str(peak_counter)+' peaks not at junctions'
    print 'Of these:'
    print str(ann_peak_counter)+' peaks at annotated sites'
    print str(GU)+' with 5prime splice sites'
    print str(AG)+' with 3prime splice sites'
    
    return rec_dict

def rec_seq(jct_dict, gff3_file, fasta_file):
    rec_dict = {}
    for transcript, junctions in jct_dict.iteritems():
        for junction in junctions:
            if junction[3] == "Recursive":
                if transcript not in rec_dict:
                    rec_dict[transcript] = [[],[]]
                rec_dict[transcript][0].append(junction[0])
                rec_dict[transcript][1].append(junction[1])
    seq_dict = SPScores.get_sequence(rec_dict, gff3_file, fasta_file)
    return seq_dict