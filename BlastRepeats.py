import pandas as pd
from subprocess import call
import shutil
import os
import sys


rgb = ["2,181,160","255,211,92","234,62,112",
       "95,135,85","1,128,181","255,130,1",
       "198,44,58","159,69,103","245,116,97",
       "75,196,213","243,182,177","171,173,89",
       "169,221,199","199,227,220","150,147,178",
       "170,66,37","140,30,43","254,197,45",
       "0,126,135","33,64,95","80,57,49",
       "180,176,173","194,187,169","127,255,212",
       "25,25,112","216,191,216","205,133,63"]

def blast_repeats(fasta, ref_fasta, out, min_e=0.01, min_coverage=0.5, top_n_results='All'):
    db_script = ['makeblastdb', '-in '+ref_fasta, '-parse_seqids', '-dbtype nucl']
    db_script = ' '.join(db_script)
    call(db_script, shell=True)

    script = ["blastn","-db "+ref_fasta,
              "-query {0}".format(fasta),'-task "blastn"',"-out {0}.out".format(out), "-outfmt 6"]
    script = ' '.join(script)
    call(script, shell=True)
    
    columns = ['Query','Contig','Identical','Length','Mismatches','Gaps','Q_start','Q_end','Ref_start','Ref_end','evalue','bitscore']
    results = pd.read_csv(out+'.out', sep='\t', names=columns)
    results.to_csv(out+'.out', sep='\t', index=False)
    
    rep_dict = {}
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                name = line.split()[0].lstrip('>').strip()
                rep_dict[name] = ''
            else:
                rep_dict[name] = rep_dict[name]+line.strip()
                
    len_dict = {k:len(v) for k,v in rep_dict.items()}
    
    for query, length in len_dict.iteritems():
        results.loc[results[results['Query'] == query].index,'Q_length'] = length
        
    results.loc[:,'Coverage'] = results['Length']/results['Q_length']
    results = results[(results['evalue'] <= min_e) & (results['Coverage'] >= min_coverage)]
    results = results.sort_values('Coverage', ascending=False).reset_index(drop=True)
    
    print "Number of hits"
    print "--------------"
    
    results.loc[:, 'Frequency'] = 0
    for name in rep_dict.keys():
        freq = len(results[results['Query'] == name])
        print name+": "+str(freq)
        results.loc[results[results['Query']==name].index, 'Frequency'] = freq
   
    if top_n_results == 'All':
        results = results.sort_values('Frequency', ascending=False)
    elif type(top_n_results) == int:
        top_results = pd.DataFrame(columns=results.columns)
        for name in set(results['Query']):
            top_n = results[results['Query'] == name].iloc[:top_n_results,:]
            top_results = top_results.append(top_n)
        results = top_results
    
    results = results.reset_index(drop=True)
    return results

def gff3_from_blast(blast_df, gff3_name, old_gff3=None):
    try:
        os.remove(gff3_name)
    except OSError:
        pass
    
    if old_gff3 is not None:
        shutil.copy(old_gff3, gff3_name)
    with open(gff3_name, 'a') as fout:
        for ix, r in blast_df.iterrows():
            if r['Ref_start'] < r['Ref_end']:
                line = [r['Contig'],'RepBase','gene',
                        str(r['Ref_start']),str(r['Ref_end']),
                        '.','+','.','ID='+r['Query']+'_'+str(ix),';\n']
            else:
                line = [r['Contig'],'RepBase','gene',
                        str(r['Ref_end']),str(r['Ref_start']),
                        '.','-','.','ID='+r['Query']+'_'+str(ix),';\n']
            line = '\t'.join(line)
            fout.write(line)
            
def bed_from_blast(blast_df, bed_name, old_gff3=None): 
    if old_gff3 is not None:
        with open(old_gff3.split('/')[-1].split('.g')[0]+'.bed', 'w') as fout:
            if old_gff3 is not None:
                fout.write('track name="Original annotation" itemRgb="On"\n')
                with open(old_gff3) as f:
                    for line in f:
                        cols = line.split('\t')
                        if len(cols) == 9:
                            if cols[2] == 'mRNA' or cols[2] == 'centromere' or cols[2] == 'telomere':
                                line = [cols[0], cols[3], cols[4], 
                                        cols[8].split('ID=')[-1].split(';')[0], '0', cols[6],
                                        cols[3], cols[4],"19,27,29",'\n']
                                line = '\t'.join(line)
                                fout.write(line)
        
    elements = list(blast_df.drop_duplicates(subset=['Query','Frequency'])['Query'])
    
    with open(bed_name, 'w') as fout:    
        fout.write('track name="New elements" itemRgb="On"\n')
        for ix, r in blast_df.iterrows():
            color = rgb[elements.index(r['Query']) % len(rgb)]
            if r['Ref_start'] < r['Ref_end']: 
                strand = '+'
                line = [r['Contig'], 
                        str(r['Ref_start']), 
                        str(r['Ref_end']), 
                        r['Query']+'_'+str(ix),
                        '0',
                        strand,
                        str(r['Ref_start']), 
                        str(r['Ref_end']), 
                        color, '\n']
            else: 
                strand = '-'
                line = [r['Contig'], 
                        str(r['Ref_end']), 
                        str(r['Ref_start']), 
                        r['Query']+'_'+str(ix),
                        '0',
                        strand,
                        str(r['Ref_end']), 
                        str(r['Ref_start']),
                        color, '\n']
            line = '\t'.join(line)
            fout.write(line)

##
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N', 'Y':'R', 'R':'Y'} 
    bases = list(seq[::-1]) 
    bases = [complement[base] for base in bases] 
    s = ''.join(bases)
    return s
        
def fasta_from_blast(blast_df, fasta_file, fasta_name):
    fa_dict = {}
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                name = line.split(' ')[0].strip('>').strip()
                fa_dict[name] = ''
            else:
                fa_dict[name] = fa_dict[name]+line.strip()
    
    with open(fasta_name, 'w') as fout:
        for ix, r in blast_df.iterrows():
            if r['Ref_start'] < r['Ref_end']:
                name = '>'+r['Query']+' '+r['Contig']+':'+str(r['Ref_start'])+'-'+str(r['Ref_end'])
                sequence = fa_dict[r['Contig']][r['Ref_start']:r['Ref_end']]
            else:
                name = '>'+r['Query']+' '+r['Contig']+':'+str(r['Ref_end'])+'-'+str(r['Ref_start'])
                sequence = fa_dict[r['Contig']][r['Ref_end']:r['Ref_start']]
                sequence = reverse_complement(sequence)
            
            fout.write(name+'\n')
            fout.write(sequence+'\n')
            
    
def main():
    if sys.argv[1] == '--help' or sys.argv[1] == '-h':
        print "Usage:"
        print "python BlastRepeats.py --query_fasta fasta_path --ref_fasta fasta_path --out output_prefix [--gff3 original_gff3] [--min_e minimum_evalue] [--min_cov minimum_coverage] [--top_results n] [--make_fasta genome_fasta_sequence]"
        print "--gff3 default: None"
        print "--min_e default: 0.01"
        print "--min_cov default: 0.5"
        print "--top_results default: 'All', provide an integer to only return the top few blast results for each query"
        print "--make_fasta, will output a fasta with an entry for each hit, must provide a fasta file for your genome."
        return None
        
    fasta = sys.argv[sys.argv.index('--query_fasta')+1]
    print fasta
    ref = sys.argv[sys.argv.index('--ref_fasta')+1]
    print ref
    out = sys.argv[sys.argv.index('--out')+1]
    print out+'\n'
    
    gff3 = None
    min_e = 0.01
    min_cov = 0.5
    fasta_file = None
    top_n_results = 'All'
    if '--gff3' in sys.argv: gff3 = sys.argv[sys.argv.index('--gff3')+1]
    if '--min_e' in sys.argv: min_e = float(sys.argv[sys.argv.index('--min_e')+1])
    if '--min_cov' in sys.argv: min_cov = float(sys.argv[sys.argv.index('--min_cov')+1])
    if '--make_fasta' in sys.argv: fasta_file = sys.argv[sys.argv.index('--make_fasta')+1]
    if '--top_results' in sys.argv: top_n_results = int(sys.argv[sys.argv.index('--top_results')+1])
        
    blast_df = blast_repeats(fasta, ref, out, min_e=min_e, min_coverage=min_cov, top_n_results=top_n_results)
    
    gff3_name = out+'.gff3'
    gff3_from_blast(blast_df, gff3_name, old_gff3=gff3)
    bed_name = out+'.bed'
    bed_from_blast(blast_df, bed_name, old_gff3=gff3)
    if fasta_file is not None:
        fasta_name = out+'.fasta'
        fasta_from_blast(blast_df, fasta_file, fasta_name)

if __name__ == "__main__":
    main()