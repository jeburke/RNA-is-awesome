import pandas as pd
import os
import numpy as np
from subprocess import Popen, PIPE
import gzip
from multiprocessing import Pool
import argparse
import pysam

def read_genome_fasta(fa):
    fa_dict = {}
    with open(fa) as f:
        for line in f:
            if line.startswith('>'):
                contig = line.split('>')[-1].strip()
                if contig not in fa_dict:
                    fa_dict[contig] = ''
                else:
                    raise KeyError("Redundant reference names in genome fasta: "+contig)
            else:
                fa_dict[contig] = fa_dict[contig]+line.strip()
    
    #for k,v in fa_dict.items():
    #    print k+': '+str(len(v))+' bp'
    
    return fa_dict

# Create index file
def build_custom_index(file_name, file_format='fasta', genome_fasta=None, flank_size=200):
    base = file_name.split('/')[-1].split('.')[0]
    if 'index' not in os.listdir('./'):
        os.makedirs('index')
    
    # First option: provide fasta file (can be compressed)
    if file_format == 'fasta':
        if not file_name.endswith('gz'):
            args = ["bowtie-build",file_name,'index/'+base]
            p = Popen(args)
        else:
            with open('fa.tmp','w') as fout:
                with gzip.open(file_name, 'rb') as f:
                    file_content = f.read()
                    for line in file_content:
                        fout.write(line)
            args = ["bowtie-build","fa.tmp",'index/'+base]
            p = Popen(args)

    # Second option: provide a table of genome coordinates
    # Must have columns Chromosome, 5p boundary, 3p boundary, index should be CNAG IDs
    else:
        if file_format == 'tsv': df = pd.read_table(file_name, index_col=0)
        if file_format == 'csv': df = pd.read_csv(file_name, index_col=0)
        else:
            raise ValueError("Unrecognized file_format: "+file_format)

        fa = read_genome_fasta(genome_fasta)
            
        with open('fa.tmp','w') as fout:
            for ix, r in df.iterrows():
                chrom, start, end = r[['Chromosome','5p boundary','3p boundary']]
                start, end = (int(start),int(end))
                lseq = fa[chrom][start-flank_size:start]
                rseq = fa[chrom][end:end+flank_size]
                
                fout.write('>'+ix+'L\n')
                fout.write(lseq+'\n')
                fout.write('>'+ix+'R\n')
                fout.write(rseq+'\n')
        
        args = ["bowtie-build","fa.tmp",'index/'+base]
        p = Popen(args)
    p.wait()
    if 'fa.tmp' in os.listdir('./'):
        os.remove('fa.tmp')
    return './index/'+base

def add_UMI((UMI_fq, read1_fq, read2_fq)):
    new_read1 = read1_fq.split('.f')[0]+'_umi.fq.gz'
    new_read2 = read2_fq.split('.f')[0]+'_umi.fq.gz'
    
    # Open new fastq for writing in zipped format
    with gzip.open(new_read1, 'wb') as fout:
        # Open UMI file
        with gzip.open(UMI_fq, 'rb') as umi:
            # Open read 1 file
            with gzip.open(read1_fq, 'rb') as f2:
                # Iterate through UMI and read 1 files together and count lines
                for n, (u, r1) in enumerate(zip(umi, f2)):
                    # 1st out of every 4 lines has read name - remember it from read 1
                    if n%4 == 0:
                        new_name = r1.strip()
                    # 2nd out of ever 4 lines has read sequence - grab from the UMI and add to read 1 name
                    elif n%4 == 1:
                        new_name = new_name+'_'+u
                        fout.write(new_name)
                        fout.write(r1)
                    else:
                        fout.write(r1)
                        
    with gzip.open(new_read2, 'wb') as fout:
        # Open UMI file
        with gzip.open(UMI_fq, 'rb') as umi:
            # Open read 1 file
            with gzip.open(read2_fq, 'rb') as f2:
                # Iterate through UMI and read 1 files together and count lines
                for n, (u, r2) in enumerate(zip(umi, f2)):
                    # 1st out of every 4 lines has read name - remember it from read 2
                    if n%4 == 0:
                        new_name = r2.strip()
                    # 2nd out of ever 4 lines has read sequence - grab from the UMI and add to read 2 name
                    elif n%4 == 1:
                        new_name = new_name+'_'+u
                        fout.write(new_name)
                        fout.write(r2)
                    else:
                        fout.write(r2)

    return (new_read1, new_read2)

def trim_filter((read1_fq, read2_fq), bound_seq1="GCCGCATCCCTGCATCCAAC", bound_seq2="CGCCTAGCAGCGGATCCAAC"):
    out1 = read1_fq.split('.f')[0]+'_trim.fq'
    out2 = read2_fq.split('.f')[0]+'_trim.fq'
    
    # Trim for right flank
    args1 = "cutadapt -g A -G {0} --trimmed-only -m 15 -o {1} -p {2} {3} {4}".format(bound_seq1, out1, out2, read1_fq, read2_fq)
    #args1 = "cutadapt -G {0} --trimmed-only -m 15 -o {1} -p {2} {3} {4}".format(bound_seq1, out1, out2, read1_fq, read2_fq)
    p = Popen(args1.split(' '))
    p.wait()
    
    # Trim for left flank
    args1 = "cutadapt -g A -G {0} --trimmed-only -m 15 -o {1} -p {2} {3} {4}".format(bound_seq2, out1+'B', out2+'B', read1_fq, read2_fq)
    #args1 = "cutadapt -G {0} --trimmed-only -m 15 -o {1} -p {2} {3} {4}".format(bound_seq2, out1+'B', out2+'B', read1_fq, read2_fq)
    p = Popen(args1.split(' '))
    p.wait()
    
    # Combine trimmed output
    for f in (out1, out2):
        with open(f,'a') as f1:
            with open(f+'B') as f2:
                for line in f2:
                    f1.write(line)
        os.remove(f+'B')
    
    return (out1, out2)

def call_bowtie(mate, index, threads=1):
    # Output name for bowtie
    sam_name = mate[0].split('_umi_trim')[0]+'.sam'
    
    # Arguments for callign bowtie
    args = ["bowtie", "-v2", "-m1", "--fr", "-X500", "-p"+str(threads), index, "-1", mate[0], "-2", mate[1],
            "--sam", sam_name]
    #print ' '.join(args)
        
    p = Popen(args, stdout=PIPE, stderr=PIPE)
    p.wait()
    for line in p.stdout:
        print(line.strip())
    for line in p.stderr:
        print(line.strip())
    return sam_name

def collapse_and_count_reads(sam):
    ref_sets = {}
    reads = pysam.Samfile(sam)
    for r in reads:
        if r.is_read1:
            if r.reference_name not in ref_sets:
                ref_sets[r.reference_name] = set()
            UMI = r.query_name.split('_')[-1]
            read_info = (r.reference_start, r.reference_end, UMI)
            ref_sets[r.reference_name].add(read_info)
            
    ref_counts = pd.DataFrame.from_dict({k:len(v) for k,v in ref_sets.items()}, orient='index')
    # Need to name column that's generated
    ref_counts = ref_counts.rename(columns={0:sam.split('_R1')[0].split('/')[-1]})
    return ref_counts

def remove_files(directory):
    for f in os.listdir(directory):
        if f.endswith('_umi.fq.gz') or f.endswith('_umi.fastq.gz'): os.remove(directory+f)
        if f.endswith('_umi_trim.fq') or f.endswith('_umi_trim.fq'): os.remove(directory+f)
        if f.endswith('.sam'): os.remove(directory+f)
    return "Finished"

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help='sub-command help', dest='subcommand')

    parser_index = subparsers.add_parser('build_index', help='Create Bowtie index')
    parser_index.add_argument("--input", help="Sequence file for creating Bowtie index")
    parser_index.add_argument("--format", default='fasta', help="Format of sequence file")
    parser_index.add_argument("--genome_fasta", default=None, help="Whole genome fasta for building index from coordinate table")
    parser_index.add_argument("--flank_size", type=int, default=200, help="Size of flanking sequences to use for index")
    
    parser_align = subparsers.add_parser('align', help='Align and count reads')
    parser_align.add_argument("--index", help="Bowtie index prefix (from build_index option)")
    parser_align.add_argument("--directory", default='./', help="Working directory containing fastq files")
    parser_align.add_argument("--threads", type=int, default=1, help="Number of processors")
    parser_align.add_argument("--out", type=str, default='read_counts', help="Name for output file")
    
    args = parser.parse_args()
    
    if args.subcommand == 'build_index':
        index_prefix = build_custom_index(args.input, file_format=args.format, 
                                          genome_fasta=args.genome_fasta, flank_size=args.flank_size)
        print("Index created at: "+index_prefix)

    elif args.subcommand == 'align':
        #Find all file sets in the directory
        #### May need to change how this done
        fq_sets = []
        for file in os.listdir(args.directory):
            if file.endswith('fq.gz') or file.endswith('fastq.gz'):
                if '_R1' in file:
                    r1 = file
                    r2 = r1.split('_R1')[0]+'_R2'+r1.split('_R1')[1]
                    ### May need to change depending on how index fastq is named
                    umi = r1.split('_R1')[0]+'_I2'+r1.split('_R1')[1]
                    fq_sets.append((args.directory+umi,args.directory+r1,args.directory+r2)) 
        
        # Add UMI
        print("Adding UMI to read names...\n")
        p = Pool(args.threads)
        barcoded = p.map(add_UMI, fq_sets)
        
        # Trim with cutadapt
        print("Trimming and filtering with Cutadapt...\n")
        trimmed = p.map(trim_filter, barcoded)
        
        # Run bowtie
        print("Aligning with Bowtie...")
        sams = []
        for pair in trimmed:
            sam = call_bowtie(pair, args.index, threads=args.threads)
            sams.append(sam)
        
        # Collapse and filter reads
        print("\nCollapsing duplicate reads and counting reads per KO flank...")
        
        all_counts = p.map(collapse_and_count_reads, sams)
        
        print('\nReads after collapsing')
        print('----------------------')
        count_df = None
        for df in all_counts:
            print(df.columns[0]+': '+str(sum(df.iloc[:,0])))
            if count_df is None:
                count_df = df
            else:
                count_df = count_df.join(df, how='outer')
        count_df = count_df.replace(np.NaN, 0)
        count_df.to_csv(args.out+'.csv')
        

        remove_files(args.directory)

if __name__ == '__main__':
    main()