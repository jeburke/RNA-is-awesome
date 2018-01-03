import sys
import pandas as pd
import os

def read_bedgraph(bedgraph):
    df = pd.read_csv(bedgraph, sep='\t', header=None)
    return df    

def write_bedgraph(dataframe, name):
    dataframe.to_csv(name+'combined.bedgraph', index=False, header=False, sep='\t')

def main(directory=None):
    bg_pairs = []
    if directory is None:
	directory = sys.argv[1]
        if len(sys.argv) == 2:
	    for file in os.listdir(directory):
	        if file.endswith('plus.bedgraph'):
#		    print file
		    bg_pairs.append((file,file.split('plus.bedgraph')[0]+'minus.bedgraph'))
        else:
	    bg_pairs.append((sys.argv[2], sys.argv[3]))
    else:
        for file in os.listdir(directory):
            if file.endswith('plus.bedgraph'):
#                print file
                bg_pairs.append((file,file.split('plus.bedgraph')[0]+'minus.bedgraph'))

    for pair in bg_pairs:
        name = pair[0].split('.bedgraph')[0].split('plus')[0]
        if not name.endswith('_'):
	    name = name+'_'
        plus = read_bedgraph(pair[0])
        minus = read_bedgraph(pair[1])
        minus[3] = minus[3].multiply(-1)

        new = plus.append(minus)
        new = new.sort_values([0,1])
        write_bedgraph(new, name)

if __name__ == "__main__":
    main() 
