'''
Usage:
    get_CNAG.py [--transcript=None] [--transcript_file=None] [--flanks]
    
Example with single transcript including 1kb upstream and downstream of TSS and PAS:
    get_CNAG.py --transcript="CNAG_01733T0" --flanks
    
Example with multiple transcripts in file with only sequence within transcript:
    get_CNAG.py --transcript_file="My_fav_CNAGs.txt"
    
    where My_fav_CNAGs.txt is a text file with a transcript ID on each new line:
        CNAG_01733T0
        CNAG_00081T0
        CNAG_04387T0
        
To run for all transcripts in the genome:
    get_CNAG.py --transcript="All"
'''

__author__ = 'jordanburke'
import re
import sys
from docopt import docopt

base_dir = sys.argv[0].split("get_CNAG")[0]
print base_dir

arguments=docopt(__doc__,version="get_CNAG 0.1")
print arguments

transcript = arguments["--transcript"]
transcript_file = arguments["--transcript_file"]

flanks = False
if arguments["--flanks"]:
    flanks = True


transcript_list = []
if transcript == "All":
    transcript_list = None
elif transcript is not None:
    transcript_list.append(transcript)

if transcript_file is not None:
    with open(transcript_file, "r") as fin:
        for line in fin:
            transcript_list.append(line.strip())
            
print transcript_list
print flanks

def fix_fasta_file(fasta_file):
    #Read fasta file for chromosome into list
        
    fasta_dict = {}
    with open(fasta_file, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                chr_num = int(line[4:6].strip())
                print chr_num
                fasta_dict[chr_num] = []
            else:
                fasta_dict[chr_num].append(line.strip())
    
    for chromosome, seq_list in fasta_dict.iteritems():
        fasta_dict[chromosome] = "".join(seq_list)
  
    with open("{0}_fixed.fa".format(fasta_file.split(".")[0]), "w") as fout:
        for chromosome, seq_list in fasta_dict.iteritems():
            fout.write(">chr"+str(chromosome)+" of Cryptococcus neoformans var. grubii H99 (CNA3)\n")
            fout.write(seq_list+"\n")
            
def load_fasta(fixed_fasta_file):
    fasta_dict = {}
    
    with open(fixed_fasta_file, "r") as fasta:
        for line in fasta:
            if line.startswith(">"):
                chr_num = int(line[4:6].strip())
            else:
                fasta_dict[chr_num]=line.strip()
    return fasta_dict
    
def load_gff3(gff3_file):
    transcript_dict = {}
    with open(gff3_file, "r") as gff3:
        for line in gff3:
            columns = re.split(r'\t+', line)
            if len(columns) > 1:
                if columns[2] == "mRNA":
                    CNAG = columns[8]
                    CNAG = CNAG[3:15]
                
                    #Transcript dictionary: keys are CNAG, values are [start, end, strand, chromosome number]
                    transcript_dict[CNAG] = [int(columns[3]), int(columns[4]), columns[6], int(columns[0][3:5].strip())]
    return transcript_dict
                
def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
        return complement(s[::-1])
    
def genbank_header(gff3_file, transcript_dict, transcript, flanks):
    CDS = []
    exons = [[],[]]
    start = transcript_dict[transcript][0]
    stop = transcript_dict[transcript][1]

    with open(gff3_file, "r") as gff3:
        for line in gff3:
            columns = re.split(r'\t+', line)
            if len(columns) > 1:
                if columns[2] == "CDS":
                    CNAG = columns[8].strip()
                    CNAG = CNAG[-12:]                 
                    if CNAG == transcript:
                        CDS = [int(columns[3]),int(columns[4])]
                
                elif columns[2] == "exon":
                    CNAG = columns[8].strip()
                    CNAG = CNAG[-12:]                   
                    if CNAG == transcript:
                        exons[0].append(int(columns[3]))
                        exons[1].append(int(columns[4]))
                    
    
    #chr1	CNA3_FINAL_CALLGENES_1	CDS	5494	5645	.	-	.	ID=CNAG_04548T0.cds;Parent=CNAG_04548T0
    #chr1	CNA3_FINAL_CALLGENES_1	exon	5322	5422	.	-	.	ID=7000009599263695;Parent=CNAG_04548T0
    
    header_list = []
    header_list.append('{:12}'.format("LOCUS")+"Chromosome "+str(transcript_dict[transcript][3])+"\t"+str(start)+"-"+str(stop))
    header_list.append('{:12}'.format("ACCESSION")+transcript)
    #header_list.append('{:12}'.format("SOURCE"))
    #header_list.append('{:12}'.format("\tORGANISM\t Cryptococcus neoformans var. grubii H99")
    header_list.append('{:12}'.format("FEATURES")+ "Location/Qualifiers")
    
    if flanks == False:
        header_list.append('{6:22}'.format("gene")+str(start-start)+".."+str(stop-start))
        header_list.append('{:23}'.format('/gene="')+transcript+'"')
        header_list.append('{6:22}'.format("CDS")+str(CDS[0]-start)+".."+str(CDS[1]-start))
        header_list.append('{:23}'.format('CDS="')+transcript+' CDS"')
                                        
    elif flanks == True:
        header_list.append('{6:22}'.format("gene")+str(start-start+1500)+".."+str(stop-start+1500))
        header_list.append('{:23}'.format('/gene="')+transcript+'"')
        header_list.append('{6:22}'.format("CDS")+str(CDS[0]-start+1500)+".."+str(CDS[1]-start+1500))
        header_list.append('{:23}'.format('CDS="')+transcript+' CDS"')
                                        
    header_list.append("ORIGIN\n")
    return "\n".join(header_list)

def get_transcript_sequence(fixed_fasta_file, gff3_file, transcript_list=None, flanks=False):
    fasta_dict = load_fasta(fixed_fasta_file)
    transcript_dict = load_gff3(gff3_file)
    
    seq_dict = {}
    
    if transcript_list is not None:
        for transcript in transcript_list:
            if len(transcript) != 12:
                print "Transcript name must be given in this format: CNAG_00001T0"
                sys.exit()
                
            start = transcript_dict[transcript][0]-1
            end = transcript_dict[transcript][1]
            strand = transcript_dict[transcript][2]
            chromosome = transcript_dict[transcript][3]
            if flanks == True:
                start = start-1500
                end = end+1500
            
            seq_dict[transcript] = fasta_dict[chromosome][start:end]
            
            if strand == "-":
                seq_dict[transcript] = reverse_complement(seq_dict[transcript])
            
            with open("{0}.fasta".format(transcript), "w") as fout:
#                fout.write(genbank_header(gff3_file, transcript_dict, transcript, flanks))
                fout.write(">"+transcript+"\n")
                fout.write(seq_dict[transcript])
#                counter = 1
                
#                while counter < len(seq_dict[transcript]):
#                    if len(seq_dict[transcript])-counter >= 60:
                        
#                        line = ('{:>9}'.format(counter)+" "+seq_dict[transcript][counter:counter+10]+" "+seq_dict[transcript][counter+11:counter+20]+" "+seq_dict[transcript][counter+21:counter+30]+" "+seq_dict[transcript][counter+31:counter+40]+" "+seq_dict[transcript][counter+41:counter+50]+" "+seq_dict[transcript][counter+51:counter+60]+"\n")
#                        fout.write(line)
#                        counter += 60
#                    else:
#                        remaining_seq = len(seq_dict[transcript])-counter
#                        num_blocks = remaining_seq/10
#                        remainder = remaining_seq%10
#                        n = 1
#                        line_list = ['{:>9}'.format(counter)]
#                        while n <= num_blocks:   
#                            line_list.append(seq_dict[transcript][counter:counter+10])
#                            n += 1
#                            counter += 10
#                        line_list.append(seq_dict[transcript][counter:counter+remainder])
#                        counter += remainder
#                        line = " ".join(line_list)
#                        fout.write(line+"\n")
#                fout.write("//")
                
    
    else:
        for transcript, coords in transcript_dict.iteritems():
            start = transcript_dict[transcript][0]
            end = transcript_dict[transcript][1]
            strand = transcript_dict[transcript][2]
            chromosome = transcript_dict[transcript][3]
            
            if flanks == True:
                start = start-1500
                end = end+1500
            
            seq_dict[transcript] = fasta_dict[chromosome][start:end]
            
            if strand == "-":
                seq_dict[transcript] = reverse_complement(seq_dict[transcript])
            
            with open("{0}.fa".format(transcript), "w") as fout:
                fout.write(seq_dict[transcript])
                
get_transcript_sequence(base_dir+"CNA3-gobs_fixed.fa", base_dir+"CNA3_FINAL_CALLGENES_1_gobs.gff3", transcript_list, flanks)       
    