import sys
sys.path.append('/home/jordan/CodeBase/RNA-is-awesome/')
sys.path.append('/home/jordan/RNA-is-awesome/')
import SPTools as SP
import json
import numpy as np
import pandas as pd
from subprocess import call
import os

def main():
    '''Usage: run SP_pipeline.py config_file untagged_sample_name organism'''
    junc_beds = []
    branch_bams = []
    CP_out = []
    CP_untagged = None
    
    with open(sys.argv[1], 'r') as config:
        for line in config:
            if 'junctions.bed' in line.lower():
                junc_beds.append(line.strip())
            elif 'branch' in line.lower():
                branch_bams.append(line.strip())
            elif sys.argv[2] in line:
                CP_untagged = line.strip()
            elif 'changepoint' in line.lower() or 'peak' in line.lower():
                CP_out.append(line.strip())

    name = sys.argv[1].split('/')[-1].split('_config.')[0]
    base_dir = sys.argv[1].split(name)[0]
    if base_dir == '': base_dir = './'
    print "Output file location and prefix: "+base_dir+name
    
    print "\nJunction bed files"
    print junc_beds
    print "\nBranch bam files"
    
    if len(branch_bams) == 2:
        print branch_bams
        use_branches = True
    elif len(branch_bams) == 0:
        print "No data for branches, continuing with only junctions"
        use_branches = False
    
    print "\nUntagged peaks"
    print CP_untagged
    print "\nChangepoint peaks"
    print CP_out
    print ''
    
    if CP_untagged is None:
        print "\n Error: no untagged file indicated"
        return None
    
    organism = sys.argv[3]
    organism, gff3, fa_dict, bowtie_index = SP.find_organism_files(organism)

    # Load in junctions
    junc_df1 = SP.build_junction_df(junc_beds[0], gff3, fa_dict, organism=organism)
    junc_df2 = SP.build_junction_df(junc_beds[1], gff3, fa_dict, organism=organism)
    
    junc_df = SP.combine_junctions(junc_df1, junc_df2)
    #print junc_df
    
    # Load in peaks
    peak_df = SP.peak_to_seq_pipeline(CP_untagged, CP_out[0], CP_out[1], gff3, fa_dict, name=name+'_CP_peaks')
    #print peak_df

    # Compare peaks and junctions
    peaks_w_junc = SP.compare_peak_junc_df(peak_df, junc_df, organism=organism)
    peaks_w_junc = SP.score_peaks(peaks_w_junc, gff3, fa_dict)
    
    # Reformat dataframe
    peaks_w_junc.index = peaks_w_junc['genome coord']
    peaks_w_junc['type index'] = np.where(peaks_w_junc['junction type'] == 'Annotated', 0, 1)
    peaks_w_junc = peaks_w_junc.sort_values('type index')
    peaks_w_junc.groupby(peaks_w_junc.index).first()
    peaks_w_junc = peaks_w_junc.drop(['index', 'type index'], axis=1)
    
    print "\nPeaks with corresponding exon-exon junctions:"
    print len(peaks_w_junc)
    print str(len(set(peaks_w_junc[~peaks_w_junc['type'].str.contains('prime')]['genome coord'])))+" unpredicted"
    
    peaks_w_junc.to_csv(base_dir+name+'_peaks_w_junc.csv')
    peaks_w_junc.to_pickle(base_dir+name+'_peaks_w_junc.pickle')
    
    # Load in branches - run main program in SPBranches to get the right output
    if use_branches is True:
        branches = SP.list_branch_points(branch_bams[0], gff3, fa_dict, organism=organism)
        branch_df = SP.create_branch_df(branches, gff3, fa_dict, organism=organism)
        if len(branch_bams) == 2:
            branches2 = SP.list_branch_points(branch_bams[1], gff3, fa_dict, organism=organism)
            branch_df2 = SP.create_branch_df(branches2, gff3, fa_dict, organism=organism)
            branch_df = branch_df.append(branch_df2)

            bed1 = branch_bams[0].split('_sorted.bam')[0]+'.bed'
            bed2 = branch_bams[1].split('_sorted.bam')[0]+'.bed'
            cat_args = "cat {0} {1} > {2}_all_branches.bed".format(bed1, bed2, name)
            call(cat_args, shell=True)

            os.remove(bed1)
            os.remove(bed2)
    
        # Compare peaks and branches
        peaks_w_branch = branch_df[branch_df['genome coord'].isin(peak_df['genome coord'])]
        peaks_w_branch = peaks_w_branch.merge(peak_df[['type','genome coord']], right_on='genome coord', left_on='genome coord', how='left')
        peaks_w_branch.index = peaks_w_branch['branch coord']

        print "\nPeaks with corresponding branches:"
        print len(peaks_w_branch)
        print str(len(set(peaks_w_branch['genome coord'])))+" unpredicted"

        peaks_w_branch.to_csv(base_dir+name+'_peaks_w_branch.csv')
        peaks_w_branch.to_pickle(base_dir+name+'_peaks_w_branch.pickle')
    
    print "\n****Finished****"

if __name__ == "__main__":
    main()