import json
import primer3
import csv
from subprocess import call
import pandas as pd
import sys

def prep_fasta(gff3, fa_dict, flank=None, gene_list=None):
    # Find all ORFs
    CDS_dict = {}
    with open(gff3) as f:
        for line in f:
            info = line.split('\t')
            if len(info) > 2:
                if info[2] == 'CDS':
                    name = info[8].split('ID=')[1].split('.cds')[0]
                    if 'Parent' in name:
                        name = name.split('Parent=')[-1].strip()
                    if name not in CDS_dict:
                        CDS_dict[name] = [info[0], [int(info[3])], [int(info[4])], info[6]]
                    CDS_dict[name][1].append(int(info[3]))
                    CDS_dict[name][2].append(int(info[4]))
    gene_dict = {}
    for name, CDS_info in CDS_dict.iteritems():
        gene = name[:-2]
        if gene not in gene_dict:
            gene_dict[gene] = [CDS_info[0], min(CDS_info[1]), max(CDS_info[2]), CDS_info[3]]
    
    if gene_list is not None:
        gene_dict = {k:v for k,v in gene_dict.items() if k in gene_list}
    
    # Retrieve sequences
    seq_dict = {}
    for gene, info in gene_dict.iteritems():
        try:
            if flank is None:
                seq = fa_dict[info[0]][info[1]:info[2]]
            else:
                seq = fa_dict[info[0]][info[1]:info[2]+flank]
            seq_dict[gene] = seq
        except KeyError:
            #print info[0]+' not in fasta file'
            pass
    return seq_dict

def run_primer3(seq):
    if len(seq) < 300:
        return None
    elif len(seq) < 500:
        size_range = [len(seq)-120,len(seq)-20]
        #print size_range
    else:
        size_range = [450,550]
        
    out = primer3.bindings.designPrimers(
    {
        'SEQUENCE_ID': 'MH1000',
        'SEQUENCE_TEMPLATE': str(seq),
        'SEQUENCE_INCLUDED_REGION': [20,len(seq)-20]
    },
    {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [size_range],
    })
    
    # Write primer3 to text file that blast can read
    keys = []
    for n in range(5):
        keys.append(('PRIMER_LEFT_'+str(n)+'_SEQUENCE','PRIMER_RIGHT_'+str(n)+'_SEQUENCE'))
    
    with open('primer3_out.txt', 'w') as fout:
        for n, (left, right) in enumerate(keys):
            fout.write('>'+str(n)+' left\n')
            fout.write(out[left]+'\n')
            fout.write('>'+str(n)+' right\n')
            fout.write(out[right]+'\n')
    return out

def blast_primers():
    script = ["blastn","-db /home/jordan/GENOMES/CNA3-gobs.fa",
              "-query primer3_out.txt",'-task "blastn-short"',"-out results.out"]
    script = ' '.join(script)
    call(script, shell=True)
    
    bad_primers = []
    with open('results.out') as f:
        min_e = []
        for line in f:
            if line.startswith("Query="):
                primer = line.split("Query= ")[1].strip()
                min_e = []
            if line.startswith("chr"):
                try:
                    e = line.split(' ')[70].strip()
                    if e == '':
                        e = line.split(' ')[69].strip()   
                except IndexError:
                    e = line.split(' ')[69].strip() 
                if 'e' in e:
                        min_e.append(e)
            if len(min_e) > 1:
                bad_primers.append(primer)
    return bad_primers

def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N', 'Y':'R', 'R':'Y'} 
    bases = list(seq) 
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def reverse_complement(s):
        return complement(s[::-1])

def get_primers(primer3_out, seq):
    keys = None
    bad_primers = blast_primers()

    # Decide which primers to use
    if len(bad_primers) > 0:
        for n in range(5):
            if str(n)+' left' in bad_primers and str(n)+' right' in bad_primers:
                pass
            else:
                keys = ['PRIMER_LEFT_'+str(n)+'_SEQUENCE','PRIMER_RIGHT_'+str(n)+'_SEQUENCE']
                break
    else:
        keys = ['PRIMER_LEFT_0_SEQUENCE','PRIMER_RIGHT_0_SEQUENCE']

    if keys is None:
        return None, None, None
    if keys is not None:
        seq1 = primer3_out[keys[0]]
        seq2 = primer3_out[keys[1]]
        amp_size = seq.find(reverse_complement(seq2))-seq.find(seq1)
        return seq1, seq2, amp_size
    


def main():
    with open('/home/jordan/GENOMES/H99_fa.json') as f:
        fa_dict = json.load(f)

    gene_df = pd.read_csv(sys.argv[1], index_col=0)
    gene_list = gene_df.index
    
    # Generate dictionary of primers designed by primer3, also runs blast during get_primers function (see above)
    seq_dict = prep_fasta('/home/jordan/GENOMES/CNA3_all_transcripts.gff3', fa_dict)
    seq_dict = {k:v for k,v in seq_dict.items() if k in gene_list}
    missing_genes = [x for x in gene_list if x not in seq_dict and str(x).startswith('CNA')]
    primer_dict = {}
    short_genes = []
    rep_genes = []
    for gene, seq in seq_dict.iteritems():
        out = run_primer3(seq)
        if out is None:
            short_genes.append(gene)
        else:
            seq1, seq2, amp_size = get_primers(out, seq)
            if seq1 is None:
                rep_genes.append(gene)
            primer_dict[gene] = (seq1, seq2, amp_size)

    # Fill in genes that are too short for the desired amplicon size - keep track of ones that are still missing
    if len(short_genes) > 0:
        short_seq_dict = prep_fasta('/home/jordan/GENOMES/CNA3_all_transcripts.gff3', fa_dict, gene_list=short_genes, flank=500)
        for gene, seq in short_seq_dict.iteritems():
            out = run_primer3(seq)
            seq1, seq2, amp_size = get_primers(out, seq)
            primer_dict[gene] = (seq1, seq2, amp_size)
    
    if len(missing_genes) > 0:
        with open('/home/jordan/GENOMES/CNA1_fa.json') as f:
            fa_dict2 = json.load(f)
        seq_dict_old = prep_fasta('/home/jordan/GENOMES/CNA1_FINAL_CALLGENES_2.annotation_fix.gff3', fa_dict2)
        for gene in missing_genes:
            try:
                out = run_primer3(seq_dict_old[gene])
                if out is not None:
                    seq1, seq2, amp_size = get_primers(out, seq_dict_old[gene])
                    if seq1 is None:
                        print gene
                    primer_dict[gene] = (seq1, seq2, amp_size)
            except KeyError:
                print "Gene not in either annotation: "+gene
        
    print str(len(rep_genes))+' genes with repetitive sequences'
    # Generate primer ordering layout
    plate_layout = []
    for letter in ['A','B','C','D','E','F','G','H']:
        for n in range(12):
            if n < 10:
                plate_layout.append(letter+'0'+str(n+1))
                plate_layout.append(letter+'0'+str(n+1))
            else:
                plate_layout.append(letter+str(n+1))
                plate_layout.append(letter+str(n+1))

    with open('500bp_primers.csv', 'w') as fout:
        fout.write('Gene,Plate,Position,Primer\n')
        for n, gene in enumerate(gene_df.index):
            if gene == '' or str(gene) == 'nan':
                pass
            else:
                if len([x for x in gene_list if x==gene]) == 2:
                    if n == list(gene_list).index(gene):
                        index = 0
                    elif n > list(gene_list).index(gene):
                        index = 1
                    plate = gene_df.loc[gene,'DP plate'][index]
                    pos = gene_df.loc[gene,'DP position'][index]
                else:
                    plate = gene_df.loc[gene,'DP plate']
                    pos = gene_df.loc[gene,'DP position']
                try:
                    # CNAG, plate, position, primer
                    if primer_dict[gene][0] is not None:
                        line_list_f = [gene, str(plate), pos, primer_dict[gene][0], str(primer_dict[gene][2])]
                        line_list_r = [gene, str(plate), pos, primer_dict[gene][1], str(primer_dict[gene][2])]
                    else: 
                        line_list_f = [gene, str(plate), pos, 'NA', 'NA']
                        line_list_r = [gene, str(plate), pos, 'NA', 'NA']
                except KeyError:
                    line_list_f = [gene, str(plate), pos, 'NA', 'NA']
                    line_list_r = [gene, str(plate), pos, 'NA', 'NA']
                line_f = ','.join(line_list_f)
                line_r = ','.join(line_list_r)
                fout.write(line_f+'\n')
                fout.write(line_r+'\n')
            
if __name__ == "__main__":
    main()
