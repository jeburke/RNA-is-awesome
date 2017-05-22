import pandas as pd

def crypto_annotate(text_file, sep=','):
    yi = pd.read_csv('/home/jordan/GENOMES/Crypto_homologs_Spombe_S288C.csv', index_col=0)
    yi.index = [x[:-2] for x in yi.index]
    yi = yi.drop_duplicates()
    hiten = pd.read_csv('/home/chomer/Code_Base/AnnotateGenes/HitensCryptoAnnot_new.txt', sep='\t', index_col=0)
    
    df = pd.read_csv(text_file, sep=sep, index_col=0)
    if df.index[0][-2] == 'T':
        df.index = [x[:-2] for x in df.index]
    
    df = df.merge(hiten, right_index=True, left_index=True)
    df = df.merge(yi, right_index=True, left_index=True)
    
    print text_file.split('.')[-2]+'_annotated.csv'
    df.to_csv(text_file.split('.')[-2]+'_annotated.csv')
    
def pombe_annotate(text_file, sep=','):
    pombase = pd.read_csv('/home/jordan/GENOMES/POMBE/gene_association.pombase', sep='\t', skiprows=42, index_col=0, usecols=[1,2,9,10,11])
    
    df = pd.read_csv(text_file, sep=sep, index_col=0)
    df = df.merge(pombase, right_index=True, left_index=True)
    df = df.drop_duplicates()
    
    print text_file.split('.csv')[0]+'_annotated.csv'
    df.to_csv(text_file.split('.csv')[0]+'_annotated.csv')