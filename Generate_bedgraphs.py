'''Usage: python Generate_bedgraphs.py directory organism'''
import sys
import GeneTools as GT

directory = sys.argv[1]
if len(sys.argv) == 2:
    organism = 'crypto'
else:
    organism = sys.argv[2]

GT.generate_scaled_bedgraphs(directory, organism=organism)