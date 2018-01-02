'''Usage: python Generate_bedgraphs.py directory organism <--start_only> <--stranded>
Arguments in <> are optional
Organism can be crypto, pombe or cerevisiae
Include the start_only argument to map only the 5' ends of reads'''


import sys
import GeneTools as GT

directory = sys.argv[1]
organism = sys.argv[2]

start_only = False
stranded = False
for arg in sys.argv:
    if arg == "--start_only":
	start_only = True
    elif arg == "--stranded":
	stranded = True

GT.generate_scaled_bedgraphs(directory, organism=organism, start_only=start_only, stranded=stranded)
