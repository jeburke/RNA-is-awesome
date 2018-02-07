'''Usage: python Generate_bedgraphs.py directory organism <--start_only> <--stranded> <--normalize untagged_sample_name> <--smooth window>
Arguments in <> are optional
Organism can be crypto, pombe or cerevisiae
Include the start_only argument to map only the 5' ends of reads'''

import sys
import GeneTools as GT
import Combine_stranded_bedgraphs
import os

def main():
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print '\nUsage:\n python Generate_bedgraphs.py directory organism <--threads n> <--start_only> <--stranded> <--normalize untagged_sample_name> <--smooth window>'
        print 'Arguments in <> are optional'
        print "directory : ./ for current directory or specify. Sorted bam files must be in specified directory. All files wil be output to specified directory."
        print 'organism : crypto, pombe, cerevisiae or candida'
        print "--threads : number of processors to use, default=1"
        print "--start_only : Map only the 5' ends of reads"
        print "--stranded : Create separate bedgraphs for the Watson and Crick strands.\n   Disclaimer - this also produces a combined bedgraph where reads on the Crick strand are negative. Please note that areas with overlapping transcription will cancel out in this representation."
        print "--normalize : Normalize to an untagged or whole cell extract sample (untagged_sample_name)"
        print "--smooth : Smooth data using rolling mean. Provide window in base pairs (e.g. 100)\n"
        return None
    
    directory = sys.argv[1]
    if not directory.endswith('/'):
        directory = directory+'/'
    organism = sys.argv[2]
    if 'crypto' not in organism.lower() and 'pombe' not in organism.lower() and 'candida' not in organism.lower() and 'cerev' not in organism.lower():
        print "Unrecognized organism"
        return None 

    threads=1
    start_only = False
    stranded = False
    normalize = False
    smooth = False
    untagged = None
    window = None
    for n, arg in enumerate(sys.argv):
        if arg == "--threads":
            try:
                threads = int(sys.argv[n+1])
            except IndexError:
                "Must provide number of threads"
                return None
            except ValueError:
                "Threads must be an integer (or number of threads not provided)"
                return None
        elif arg == "--start_only":
            start_only = True
        elif arg == "--stranded":
            stranded = True
        elif arg == "--normalize":
            normalize = True
            try:
                untagged = sys.argv[n+1].split('/')[-1]
            except IndexError:
                print "Must provide untagged sample name"
                return None
            
            if untagged.startswith('--'):
                print "Must provide untagged sample name"
                return None
        elif arg == "--smooth":
            smooth = True
            try:
                window = int(sys.argv[n+1])
            except IndexError:
                print "Must provide window size"
                return None
            except ValueError:
                print "Must provide window size"
                return None

    print "Generating scaled bedgraphs..."
    GT.generate_scaled_bedgraphs2(directory, organism=organism, start_only=start_only, stranded=stranded, threads=threads)
    
    if stranded is True:
        Combine_stranded_bedgraphs.main(directory=directory)

    if normalize is True:
        print "\nNormalizing to untagged..."
        if untagged.endswith('.bam'):
            untagged = untagged.split('.bam')[0]
        bg_list = [directory+x for x in os.listdir(directory) if x.endswith('.bedgraph')]
        untagged_bg = [x for x in bg_list if untagged in x][0]
        bg_list.remove(untagged_bg)

        for bg in bg_list:
            GT.normalize_bedgraph(bg, untagged_bg)

    if smooth is True:
        print "\nSmoothing with {0} bp window...".format(str(window))
        bg_list = [directory+x for x in os.listdir(directory) if x.endswith('.bedgraph')]
        GT.smooth_bedgraphs(bg_list, window)
        
if __name__ == "__main__":
    main()