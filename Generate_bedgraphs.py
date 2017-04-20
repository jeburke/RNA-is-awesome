'''Usage: python Generate_bedgraphs.py directory <organism> <start_only>
Arguments in <> are optional
Organism can be crypto, pombe or cerevisiae
Include the start_only argument to map only the 5' ends of reads'''


import sys
import GeneTools as GT

directory = sys.argv[1]
if len(sys.argv) == 2:
    organism = 'crypto'
if len(sys.argv) == 3 and sys.argv[-1] == 'start_only':
    organism = 'crypto'
else:
    organism = sys.argv[2]

if sys.argv[-1] == 'start_only':
    start_only = True
else: start_only = False
    
    
GT.generate_scaled_bedgraphs(directory, organism=organism, start_only=start_only)