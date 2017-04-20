import sys
import os
import subprocess

def align_fastq(directory, threads=1, organism='crypto'):
    if directory[-1] != '/':
        directory = directory+'/'
        
    if 'crypto' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/Crypto_for_gobs'
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/S288C/S288C'
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
    elif 'pombe' in organism.lower():
        bowtie_ix = '/home/jordan/GENOMES/POMBE/Spombe'
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
        
    fastq_list = []
    for file in os.listdir(directory):
        if file.lower().endswith("fastq.gz"):
            fastq_list.append(directory+file)
            
    for fastq in fastq_list:
        prefix = fastq.split('.fastq')[0]
        args = 'tophat2 --read-mismatches 3 --read-gap-length 2 --read-edit-dist 3 --min-anchor-length 8 --splice-mismatches 1 --min-intron-length 50 --max-intron-length 1200 --max-insertion-length 3 --max-deletion-length 3 --num-threads {0} --max-multihits 100 --library-type fr-firststrand --segment-mismatches 3 --no-coverage-search --segment-length 20 --min-coverage-intron 10 --max-coverage-intron 100000 --min-segment-intron 50 --max-segment-intron 500000 --b2-sensitive --bowtie1 -G {1} -o {2} {3} {4}'.format(str(threads), gff3, prefix, bowtie_ix, fastq)
        #print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        p.wait()
        
        print '\n'
        print prefix
        with open(prefix+'/align_summary.txt','r') as f:
            for line in f:
                print line
            print '\n'
            
    print 'Finished'

def sort_index_bam_files(base_dir):
    try:
        check_samtools = subprocess.check_output(['samtools','--version'])
        version = check_samtools.split('\n')[0].split(' ')[1]
    except subprocess.CalledProcessError:
        version = 'old'
    print version
    
    
    if base_dir[-1] != '/':
        base_dir = base_dir+'/'
        
    tophat_out = []
    names = []
    for name in os.listdir(base_dir):
        if os.path.isdir(base_dir+name):
            tophat_out.append(base_dir+name+'/accepted_hits.bam')
            names.append(base_dir+name)
            print base_dir+name+'/accepted_hits.bam'
    n=0
    sorted_bams = []
    for n in range(len(tophat_out)):
        name = names[n]
        args = 'samtools sort {0} -o {1}'.format(tophat_out[n], name+'_sorted.bam')
        if version == 'old':
            args = 'samtools sort {0} {1}'.format(tophat_out[n], name+'_sorted.bam')
        print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
        sorted_bams.append(name+'_sorted.bam')
    p.wait()
        
    for bam in sorted_bams:
        args = 'samtools index {0}'.format(bam)
        print args
        p = subprocess.Popen(args.split(' '), stdout=subprocess.PIPE)
    p.wait()
        
def count_with_HTseq(base_dir, organism='crypto'):
    if base_dir[-1] != '/':
        base_dir = base_dir+'/'
    
    if 'crypto' in organism.lower():
        gff3 = '/home/jordan/GENOMES/CNA3_all_transcripts.gff3'
    elif 'cerev' in organism.lower():
        gff3 = '/home/jordan/GENOMES/S288C/saccharomyces_cerevisiae_R64-2-1_20150113.gff3'
    elif 'pombe' in organism.lower():
        gff3 = '/home/jordan/GENOMES/POMBE/schizosaccharomyces_pombe.chr.gff3'
    
    
    bam_files = []
    for file in os.listdir(base_dir):
        if file.endswith('_sorted.bam'):
            bam_files.append(file)
    
    for bam in bam_files:
        args = 'htseq-count -f bam -s reverse -t gene -i ID {0} {1}'.format(base_dir+bam, gff3)
        print args
        
        out = subprocess.check_output(args.split(' '))
        with open(base_dir+bam.split('_sorted')[0]+'.htseq','w') as fout:
            fout.write(out)
        
def main():
    # Input will be RNAseq_tools.py directory number_of_threads_tophat organism
    # Current organism options: 'crypt', 'pombe', 'cerevisiae'
    if len(sys.argv) == 3:
        organism = 'crypto'
    else:
        organism = sys.argv[3]
        
    directory = sys.argv[1]
    threads = sys.argv[2]
    
    align_fastq(directory, threads=threads, organism=organism)
    sort_index_bam_files(directory)
    count_with_HTseq(directory, organism=organism)
    
if __name__ == "__main__":
    main()