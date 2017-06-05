import pandas as pd

def crypto_annotate(text_file, sep=',', multi_index=False):
    '''Annotates a spreadsheet or text file where the gene ID (CNAG) is the first column.
    Uses Hiten's hand annotation and S. cerevisiae and S. pombe orthologs determined by Yi.
    
    Parameters
    ----------
    text_file : str
              Your spreadsheet or text file
    sep : str, default ','
              Column separator character. Change to '\t' if tab separated.
    
    Output
    -------
    csv file : File (comma separated) with annotation appended to each row.
              '''
    
    yi = pd.read_csv('/home/jordan/GENOMES/Crypto_homologs_Spombe_S288C.csv', index_col=0)
    yi.index = [x[:-2] for x in yi.index]
    yi = yi.drop_duplicates()
    hiten = pd.read_csv('/home/chomer/Code_Base/AnnotateGenes/HitensCryptoAnnot_new.txt', sep='\t', index_col=0)
    
    if multi_index is True:
        df = pd.read_csv(text_file, sep=sep, index_col=0, header=[0,1])
    else:
        df = pd.read_csv(text_file, sep=sep, index_col=0)
        
    if len(df.index) == 0:
        print "Empty spreadsheet, cannot annotate"
        return

    if df.index[0][-2] == 'T':
        df.index = [x[:-2] for x in df.index]
        
    
    df = df.merge(hiten, right_index=True, left_index=True, how='left')
    df = df.merge(yi, right_index=True, left_index=True, how='left')
    df = df.groupby(df.index).first()
    
    print text_file.split('.')[-2]+'_annotated.csv'
    df.to_csv(text_file.split('.')[-2]+'_annotated.csv')
    
def pombe_annotate(text_file, sep=',', multi_index=False):
    '''Annotates a spreadsheet or text file where the gene ID is the first column.
    Uses Pombase annotation.
    
    Parameters
    ----------
    text_file : str
              Your spreadsheet or text file
    sep : str, default ','
              Column separator character. Change to '\t' if tab separated.
    
    Output
    -------
    csv file : File (comma separated) with annotation appended to each row.
              '''
    
    pombase = pd.read_csv('/home/jordan/GENOMES/POMBE/gene_association.pombase', sep='\t', skiprows=42, index_col=0, usecols=[1,2,9,10,11])
    
    if multi_index is True:
        df = pd.read_csv(text_file, sep=sep, index_col=0, header=[0,1])
    else:
        df = pd.read_csv(text_file, sep=sep, index_col=0)
        
    if df.index[0][-2] == '.':
        df.index = [x[:-2] for x in df.index]
    
    if df.index[0].startswith('gene'):
        df.index = [x.split('gene:')[1] for x in df.index]
        
    df = df.merge(pombase, right_index=True, left_index=True, how='left')
    df = df.drop_duplicates()
    
    print text_file.split('.csv')[0]+'_annotated.csv'
    df.to_csv(text_file.split('.csv')[0]+'_annotated.csv')