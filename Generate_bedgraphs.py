'''Usage: python Generate_bedgraphs.py directory organism <--start_only> <--stranded> <--normalize untagged_sample_name> <--smooth window>
Arguments in <> are optional
Organism can be crypto, pombe or cerevisiae
Include the start_only argument to map only the 5' ends of reads'''

import sys
import GeneTools as GT
import Combine_stranded_bedgraphs
import os

directory = sys.argv[1]
organism = sys.argv[2]

start_only = False
stranded = False
normalize = False
smooth = False
untagged = None
window = None
for n, arg in enumerate(sys.argv):
    if arg == "--start_only":
        start_only = True
    elif arg == "--stranded":
        stranded = True
    elif arg == "--normalize":
        normalize = True
        try:
            untagged = sys.argv[n+1]
        except IndexError:
            print "Must provide untagged sample name"
        if untagged.startswith('--'):
            print "Must provide untagged sample name"
    elif arg == "--smooth":
        smooth = True
        try:
            window = int(sys.argv[n+1])
        except IndexError:
            print "Must provide window size"
        except ValueError:
            print "Must provide window size"

print "Generating scaled bedgraphs..."
GT.generate_scaled_bedgraphs(directory, organism=organism, start_only=start_only, stranded=stranded)
if stranded is True:
    Combine_stranded_bedgraphs.main(directory=directory)
    
if normalize is True:
    print "Normalizing to untagged..."
    bg_list = [x for x in os.listdir(directory) if x.endswith('.bedgraph')]
    untagged_bg = [x for x in bg_list if untagged in x]
    bg_list.remove(untagged_bg)
    
    for bg in bg_list:
        GT.normalize_bedgraph(bg, untagged_bg)
        
if smooth is True:
    print "Smoothing with {0} bp window...".format(str(window))
    bg_list = [x for x in os.listdir(directory) if x.endswith('.bedgraph')]
    GT.smooth_bedgraphs(bg_list, window)