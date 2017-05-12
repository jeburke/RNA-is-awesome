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
from subprocess import call
import json
import HTSeq
import itertools
from multiprocessing import Pool
from functools import partial
import os
from itertools import izip,islice,tee
import re

### Collect all 5p splice sites of interest and get 20 nt of sequence before splice site - pass in splice site dict
def collect_intron_seq(gff3_file, fasta_file, ss_dict=None, junction_bed=None, gene_list=None, peak_df=None, organism=None):
    transcript_dict = SP.build_transcript_dict(gff3_file, organism=organism)
    if type(fasta_file) == dict:
        fasta_dict = fasta_file
    elif fasta_file.endswith('json'):
        with open(fasta_file, 'r') as f:
            fasta_dict = json.load(f)
    else:
        fasta_dict = make_fasta_dict(fasta_file)
    if ss_dict is not None:
        ss_dict=ss_dict
    elif junction_bed is not None:
        ss_dict = SP.build_junction_dict(junction_bed, gff3_file, transcript_dict, organism=organism)
    elif peak_df is not None:
        ss_dict = {}
        peak_df = peak_df[~peak_df['type'].str.contains('prime')]
        for ix, r in peak_df.iterrows():
            if r['transcript'] not in ss_dict:
                ss_dict[r['transcript']] = []
            if r['strand'] == '+':
                ss_dict[r['transcript']].append((r['position'],r['position']+50))
            elif r['strand'] == '-':
                ss_dict[r['transcript']].append((r['position'],r['position']-50))
                
    else:
        ss_dict, intron_flag = SP.list_splice_sites(gff3_file, gene_list=gene_list, organism=organism)
        ss_dict = SP.collapse_ss_dict(ss_dict)
    
    print len(ss_dict)
    
    seq_dict = {}
    for transcript, introns in ss_dict.iteritems():
        if junction_bed is None:
            if organism == 'pombe':
                transcript = transcript+'.1'
            else:
                transcript = transcript+'T0'
        introns = list(introns)
        strand = transcript_dict[transcript][2]
        chrom = transcript_dict[transcript][3]
        n = 0
        for n in range(len(introns)):
            if strand == '+':
                seq_dict[transcript+'-'+chrom+':'+str(introns[n][0]+1)] = fasta_dict[chrom][introns[n][0]+2:introns[n][0]+17]
            elif strand == '-':
                seq = fasta_dict[chrom][introns[n][0]-16:introns[n][0]-1]
                seq_dict[transcript+'-'+chrom+':'+str(introns[n][0])] = SP.reverse_complement(seq)
    return seq_dict

### Samtools view and grep for each sequence
def sub_findre(s,substring,diffnumber):
    sublen=len(substring)
    zip_gen=(izip(substring,islice(s,i,i+sublen)) for i in xrange(len(s)))
    for z in zip_gen:
        l,z=tee(z)
        if sum(1 for i,j in l if i==j)>=sublen-diffnumber:
            new=izip(*z)
            next(new)
            yield ''.join(next(new))
            
def sub_findre2(s,substring):
    if substring not in sub_findre2.cache:
        sub_findre2.cache[substring] = [
            re.compile('.'.join((substring[:x], substring[x+1:])))
            for x in range(len(substring))]
    for regex in sub_findre2.cache[substring]:
        for match in regex.findall(s):
            yield match
sub_findre2.cache = {}

def a_rich(string):
    a_count = 0
    for char in string:
        if char == 'A':
            a_count += 1
    perc_A = float(a_count)/len(string)
    if perc_A >= 0.4:
        high_a = True
    else:
        high_a = False
    return high_a
            
def process_reads((fastq_file, seq_dict), prefix=None, mismatches=0):
    if prefix is None:
        prefix = fastq_file
    fa_line_list = []
    unsplit_fa_line_list = []
    
    if '.fastq' in fastq_file or '.fq0' in fastq_file or '.fq1' in fastq_file or fastq_file.endswith('fq'):
        fastq = HTSeq.FastqReader(fastq_file, "solexa")
    elif '.fasta' in fastq_file or '.fa0' in fastq_file or '.fa1' in fastq_file or fastq_file.endswith('fa'):
        fastq = HTSeq.FastaReader(fastq_file)
    
    if '_unsplit' in fastq_file:
        mismatches = 1
    
    for read in fastq:
        new_read = None
        for intron, seq in seq_dict.iteritems():
            if seq in read.seq:
                read_seq = read.seq.split(seq)[0]
                new_length = len(read_seq)
                new_read = [read.name+':'+intron, read_seq]
                break
            elif seq not in read.seq and mismatches > 0:
                in_read = list(sub_findre2(read.seq, seq))
                if len(in_read) > 0:
                    for n, subseq in enumerate(in_read):
                        read_seq = read.seq.split(subseq)[0]
                        new_length = len(read_seq)
                        new_read = [read.name+':'+intron+'-'+str(n), read_seq]
                    break
        
        if new_read is not None and len(new_read[1]) > 15:
            fa_line_list.append(('>'+new_read[0]+'\n',new_read[1]+'\n'))
        
        if new_read is None and a_rich(read.seq) is False:
            unsplit_fa_line_list.append(('>'+read.name+'\n',read.seq+'\n'))
        
        if len(fa_line_list) == 5000:
            with open('{}_split.fa'.format(prefix), 'a') as fout:
                for pair in fa_line_list:
                    fout.write(pair[0])
                    fout.write(pair[1])
            fa_line_list = []
        
        if len(unsplit_fa_line_list) == 5000:
            with open('{}_unsplit.fa'.format(prefix), 'a') as fout:
                for pair in unsplit_fa_line_list:
                    fout.write(pair[0])
                    fout.write(pair[1])
            unsplit_fa_line_list = []            
            
    with open('{}_split.fa'.format(prefix), 'a') as fout:
        for pair in fa_line_list:
            fout.write(pair[0])
            fout.write(pair[1])
            
    with open('{}_unsplit.fa'.format(prefix), 'a') as fout:
        for pair in unsplit_fa_line_list:
            fout.write(pair[0])
            fout.write(pair[1]) 
            
    with open('fa_list.txt','a') as f:
        f.write('{}_split.fa\n'.format(prefix))
        
    with open('unsplit_fa_list.txt','a') as f:
        f.write('{}_unsplit.fa\n'.format(prefix))
        

def find_split_reads(fastq_file, seq_dict, prefix, threads=0):
    if threads > 0:
        p = Pool(threads)
        fastq_list = split_fastq_file(fastq_file, threads)
        params = []
        for fastq in fastq_list:
            params.append((fastq, seq_dict))
        p.map(process_reads, params)
        p.close()
        p.join()
        
        fa_list = []
        with open('fa_list.txt','r') as f:
            for line in f:
                fa_list.append(line.strip())

        unsplit_fa_list = []
        with open('unsplit_fa_list.txt','r') as f:
            for line in f:
                unsplit_fa_list.append(line.strip())                

        call_list = ["cat",fastq_file+"*_split.fa",">","{0}_split.fa".format(prefix)]
        script = ' '.join(call_list)
        print script
        call(script, shell=True)
        
        call_list = ["cat",fastq_file+"*_unsplit.fa",">","{0}_unsplit.fa".format(prefix)]
        script = ' '.join(call_list)
        print script
        call(script, shell=True)
        
        for fastq in fastq_list:
            os.remove(fastq)
        try:
            for fa in fa_list:
                os.remove(fa)

            for fa in unsplit_fa_list:
                os.remove(fa)
        except OSError:
            pass
            
        os.remove('fa_list.txt')
        os.remove('unsplit_fa_list.txt')
        
    else:
        process_reads((fastq_file, seq_dict), prefix=prefix)
    
def split_fastq_file(fastq_file, threads):
    with open(fastq_file, 'r') as f:
        for i, l in enumerate(f):
            pass
    file_len = i + 1
    num_files = file_len/20000
    print num_files
    
    call(["split", "-d", "-l 20000" , "-a 5", fastq_file, fastq_file])
    
    n=0
    fastq_list = []
    for n in range(num_files):
        fastq_list.append(fastq_file+format(n, '05'))
    print fastq_list[-1]
    return fastq_list
        
    
def list_branch_points(sorted_bam_file, gff3_file, fasta_dict, organism=None):
    transcript_dict = SP.build_transcript_dict(gff3_file, organism=organism)
    if type(fasta_dict) == str:
        with open(fasta_dict, 'r') as f:
            fasta_dict = json.load(f)
    
    branch_dict = {}
    read_counter = 0
    br_counter = 0
    bam_reader = HTSeq.BAM_Reader(sorted_bam_file)
    for a in bam_reader:
        read_counter += 1
        transcript = a.read.name.split('-')[0].split(':')[-1]
        splice_site = a.read.name.split('-')[1]
        if transcript not in branch_dict:
            branch_dict[transcript] = {}
        if splice_site not in branch_dict[transcript]:
            branch_dict[transcript][splice_site] = []
        if a.iv is not None:
            strand = a.iv.strand
            read_end = a.iv.end
            if strand == '-':
                read_end = a.iv.start
            if strand == transcript_dict[transcript][2]:
                branch_dict[transcript][splice_site].append(read_end)
                br_counter += 1
                
    print "Reads analyzed: "+str(read_counter)
    print "Reads assigned as branches: "+str(br_counter)
    
    new_branch_dict = {}
    for transcript, introns in branch_dict.iteritems():
        new_branch_dict[transcript] = []
        for intron, branches in introns.iteritems():
            new_branch_list = []
            new_branch_counts = []
            for branch in branches:
                flag = False
                if len(new_branch_list) > 0:
                    for pos in range(branch-2,branch+3):
                        if pos in new_branch_list: 
                            flag = True
                            br_id = new_branch_list.index(pos)
                            new_branch_counts[br_id] += 1
                if flag == False: 
                    new_branch_list.append(branch)
                    new_branch_counts.append(1)
            if len(new_branch_list) > 0:
                new_branch_dict[transcript].append([intron, new_branch_list, new_branch_counts])
            
    with open('{0}_branches.bedgraph'.format(sorted_bam_file.split('.')[0]), 'w') as fout:
        for transcript, introns in new_branch_dict.iteritems():
            for intron in introns:
                chrom = intron[0].split(':')[0]
                n = 0
                for n in range(len(intron[1])):
                    start = intron[1][n]
                    end = intron[1][n]+2
                    value = intron[2][n]
                    line_list = [chrom, start, end, value, '\n']
                    line_list = map(str, line_list)
                    line = '\t'.join(line_list)
                    fout.write(line)
    
    with open('{0}_branches.bed'.format(sorted_bam_file.split('.')[0]), 'w') as fout:
        fout.write('track name=junctions description="TopHat junctions"\n')
        for transcript, introns in new_branch_dict.iteritems():
            strand = transcript_dict[transcript][2]
            for intron in introns:
                chrom = intron[0].split(':')[0]
                start = int(intron[0].split(':')[1])
                n=0
                for n in range(len(intron[1])):
                    end = intron[1][n]
                    value = intron[2][n]
                    size = abs(end-start)+30
                    if abs(end-start) > 2000:
                        #print intron
                        pass
                    elif abs(end-start) > 5:
                        #[seqname] [start] [end] [id] [score] [strand] [thickStart] [thickEnd] [r,g,b][block_count] [block_sizes] [block_locations]
                        read_id = intron[0]+'-'+str(n)
                        block_size = '0,'+str(size)
                        line_list = [chrom, str(start-30), str(end+30), read_id, str(value), strand, str(start-30), str(end+30), '75,196,213', '2', '30,30', block_size, '\n']
                        line = '\t'.join(line_list)
                        fout.write(line)
    
    with open('{0}_branches.json'.format(sorted_bam_file.split('.')[0]), 'w') as fout:
        json.dump(new_branch_dict, fout)
    
    return new_branch_dict
            

def find_multibranch_introns(branch_dict, df_size):
    columns = ['transcript','chromosome','5p splice site','branch site','reads']
    df = pd.DataFrame(columns=columns, index= range(df_size))
    n = 0
    for transcript, sites in branch_dict.iteritems():
        for site in sites:
            branch_list = site[1]
            if len(branch_list) > 1:
                m = 0
                for m in range(len(branch_list)):
                    row = [transcript, site[0].split(':')[0], site[0].split(':')[1], branch_list[m], site[2][m]]
                    df.iloc[n] = row
                    if n%10000 == 0:
                        print n
                    n += 1
    df = df.dropna()
    return df       


###Maybe remove this function. Takes too long.
def best_bp_seq(branch_dict, df_size, fasta_dict, transcript_dict):
    #Generate a dataframe with introns with multiple branch sites
    df = find_multibranch_introns(branch_dict, df_size)
    
    #Pick the best 'branch' (main issue being that sometimes RT skips the branch and there's a sequence shift)
    best_bp_df = pd.DataFrame(columns = df.columns)
    br_tx_list = list(set(df['transcript'].tolist()))
    for transcript in br_tx_list:
        tx_df = df[df['transcript'] == transcript]
        five_sites = list(set(tx_df['5p splice site'].tolist()))
        for site in five_sites:
            site_df = tx_df[tx_df['5p splice site'] == site]
            site_df = site_df[site_df['branch site'] != float(site)]
            site_df = site_df.reset_index()
            bp_indexes = set()
            no_set = set()
            bp_list = site_df['branch site'].tolist()
            for bp in bp_list:
                bp_dists = [x-bp for x in bp_list]
                flag = False
                for dist in bp_dists:
                    if dist < 0 and dist > -5:
                        no_set.add(bp_list.index(bp))
                    else:
                        bp_indexes.add(bp_list.index(bp))

            bp_indexes = bp_indexes.difference(no_set) 
            bp_indexes = list(bp_indexes)

            site_df = site_df[site_df.index.isin(bp_indexes)].sort_values('reads', ascending=False)
            site_df = site_df.reset_index()

            if len(site_df) > 0:
                best_bp_df = best_bp_df.append(site_df.ix[0])
                
    best_bp_df.drop('level_0', axis=1, inplace=True)
    best_bp_df = best_bp_df.reset_index()
    best_bp_df.drop('level_0', axis=1, inplace=True)
    
    #Fetch the sequences for the branch sites
    branch_seq_df = best_bp_df
    
    five_seq_list = []
    bp_seq_list = []
    for index, row in best_bp_df.iterrows():
        chrom = row['chromosome']
        five_site = int(row['5p splice site'])
        bp = int(row['branch site'])
        strand = transcript_dict[transcript][2]
        if strand == '+':
            five_seq_list.append(fasta_dict[chrom][five_site-2:five_site+6])
            bp_seq = fasta_dict[chrom][bp-2:bp+3]
            if bp_seq.endswith('A'):
                bp_seq = fasta_dict[chrom][bp-1:bp+4]
            bp_seq_list.append(bp_seq)
        elif strand == '-':
            five_seq = fasta_dict[chrom][five_site-6:five_site+2]
            five_seq_list.append(SPPeaks.reverse_complement(five_seq))
            bp_seq = fasta_dict[chrom][bp:bp+5]
            bp_seq = SPPeaks.reverse_complement(bp_seq)
            if bp_seq.endswith('A'):
                bp_seq = fasta_dict[chrom][bp-1:bp+4]
                bp_seq = SPPeaks.reverse_complement(bp_seq)
            bp_seq_list.append(bp_seq)
    branch_seq_df['5p sequence'] = pd.Series(five_seq_list, index = branch_seq_df.index)
    branch_seq_df['Branch sequence'] = pd.Series(bp_seq_list, index = branch_seq_df.index)
    branch_seq_df = branch_seq_df[branch_seq_df['reads'] > 10.0]

    print len(best_bp_df)
    print len(branch_seq_df)
    return branch_seq_df
    
def annotate_branch_df(branch_df, gff3, fasta_dict, gene_list=None, organism=None):
    ss_dict, intron_flag = SP.list_splice_sites(gff3, gene_list=gene_list, organism=organism)
    ss_dict = SP.collapse_ss_dict(ss_dict)
    
    #Add 3' splice sites and intron sizes
    intron_sizes = []
    three_sites = []
    strands = []
    for index, row in branch_df.iterrows():
        if type(row['transcript']) == float:
            print row
        introns = ss_dict[(row['transcript'][:-2])]
        flag = False
        for intron in introns:
            if int(row['5p splice site']) == intron[0]+1:
                intron_sizes.append(intron[1]-intron[0])
                three_sites.append(intron[1])
                strands.append('+')
                flag = True
                break
            elif int(row['5p splice site']) == intron[0]:
                intron_sizes.append(intron[0]-intron[1])
                three_sites.append(intron[1])
                strands.append('-')
                flag = True
                break
        if flag is False:
            print row['transcript']
            print row['5p splice site']
    print len(intron_sizes)
    print len(three_sites)

    new_branch_df = branch_df
    new_branch_df['intron size'] = pd.Series(intron_sizes, index= new_branch_df.index)
    new_branch_df['3p splice site'] = pd.Series(three_sites, index= new_branch_df.index)
    new_branch_df['strand'] = pd.Series(strands, index= new_branch_df.index)
    
    #Find sequeces for branches and spilce sites
    branch_seq_df = new_branch_df
    five_seq_list = []
    three_seq_list = []
    bp_seq_list = []
    
    for index, row in new_branch_df.iterrows():
        chrom = row['chromosome']
        five_site = int(row['5p splice site'])
        three_site = int(row['3p splice site'])
        bp = int(row['branch site'])
        strand = row['strand']
        if strand == '+':
            five_seq_list.append(fasta_dict[chrom][five_site-2:five_site+6])
            three_seq_list.append(fasta_dict[chrom][three_site-5:three_site+3])
            bp_seq = fasta_dict[chrom][bp-2:bp+3]
            if bp_seq.endswith('A'):
                bp_seq = fasta_dict[chrom][bp-1:bp+4]
            bp_seq_list.append(bp_seq)
        elif strand == '-':
            five_seq = fasta_dict[chrom][five_site-6:five_site+2]
            five_seq_list.append(SP.reverse_complement(five_seq))
            three_seq = fasta_dict[chrom][three_site-2:three_site+6]
            three_seq_list.append(SP.reverse_complement(three_seq))
            bp_seq = fasta_dict[chrom][bp:bp+5]
            bp_seq = SP.reverse_complement(bp_seq)
            if bp_seq.endswith('A'):
                bp_seq = fasta_dict[chrom][bp-1:bp+4]
                bp_seq = SP.reverse_complement(bp_seq)
            bp_seq_list.append(bp_seq)
    
    branch_seq_df['5p sequence'] = pd.Series(five_seq_list, index = branch_seq_df.index)
    branch_seq_df['3p sequence'] = pd.Series(three_seq_list, index = branch_seq_df.index)
    branch_seq_df['Branch sequence'] = pd.Series(bp_seq_list, index = branch_seq_df.index)
    branch_seq_df = branch_seq_df[branch_seq_df['reads'] > 10.0]
    
    return branch_seq_df

def write_seq_list_to_file(df, prefix):
    with open('{0}_5pseq_list.txt'.format(prefix),'w') as f:
        for seq in df['5p sequence'].tolist():
            f.write(seq+'\n')
    with open('{0}_3pseq_list.txt'.format(prefix),'w') as f:
        for seq in df['3p sequence'].tolist():
            f.write(seq+'\n')
    with open('{0}_BPseq_list.txt'.format(prefix),'w') as f:
        for seq in df['Branch sequence'].tolist():
            f.write(seq+'\n')
            
def find_best_branch(branch_seq_df):
    best_bp_df = pd.DataFrame(columns = branch_seq_df.columns)
    br_tx_list = list(set(branch_seq_df['transcript']))
    
    for transcript in br_tx_list:
        tx_df = branch_seq_df[branch_seq_df['transcript'] == transcript]
        five_sites = list(set(tx_df['5p splice site'].tolist()))
        for site in five_sites:
            site_df = tx_df[tx_df['5p splice site'] == site]
            site_df = site_df[site_df['branch site'] != float(site)]
            site_df = site_df.reset_index()
            bp_indexes = set()
            no_set = set()
            bp_list = site_df['branch site'].tolist()
            for bp in bp_list:
                bp_dists = [x-bp for x in bp_list]
                flag = False
                for dist in bp_dists:
                    if dist < 0 and dist > -5:
                        no_set.add(bp_list.index(bp))
                    else:
                        bp_indexes.add(bp_list.index(bp))

            bp_indexes = bp_indexes.difference(no_set) 
            bp_indexes = list(bp_indexes)

            site_df = site_df[site_df.index.isin(bp_indexes)].sort_values('reads', ascending=False)
            site_df = site_df.reset_index()

            if len(site_df) > 0:
                best_bp_df = best_bp_df.append(site_df.ix[0])

    best_bp_df.drop('level_0', axis=1, inplace=True)
    best_bp_df = best_bp_df.reset_index()
    best_bp_df.drop('level_0', axis=1, inplace=True)
    print len(best_bp_df)
    return best_bp_df