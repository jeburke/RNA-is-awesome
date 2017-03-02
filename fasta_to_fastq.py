import sys
output = sys.argv[2]
with open(output, 'w') as fout:
    with open(sys.argv[1], 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = '@'+line.strip()[1:]
            elif 'TCTTCTGCTTG' in line:
	        pass
	    else:
	        seq = line.strip()
		score = 'I' * len(seq)
		fout.write(name+'\n')
		fout.write(seq+'\n')
		fout.write('+\n')
		fout.write(score+'\n')
		
            
