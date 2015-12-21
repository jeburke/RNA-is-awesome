__author__ = 'jordanburke'

''' This script compares splice site scores and normalized tables from SPNormalizeAndPlot and plots the results
Usage: NormalizeAndPlot.py <splice site scores> <normalized exon table> <normalized intron table> <prefix> '''

import sys
import SPTools
import pandas as pd
import matplotlib.pyplot as plt

def merge_tables(n):
    file = sys.argv[n]
    counts_df = SPTools.build_intron_tables(file)
    merged_df = pd.merge(scores_df, counts_df, right_index=True, left_index=True)
    scores_df_sorted = merged_df.sort()
    fout = open("{0}_scored_ss{1}.txt".format(sys.argv[-1].split(".")[0], str(n)), "w")
    print fout
    scores_df_sorted.to_csv(fout, sep='\t')
    fout.close()
    return scores_df_sorted

def get_points(normalizedResults, a):
    values = []
    n = 0
    for name in normalizedResults.columns:
        if name.split("_")[-1] == "norm":
            if n == a:
                for entry in normalizedResults[name]:
                    values.append(entry)
                n += 1
            else:
                n += 1
    return values

def scatter_plot(x,y1,y2,label1, label2, title):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y1, c='royalblue', label=label1, alpha = 0.5, edgecolor='darkslateblue')
    ax1.scatter(x, y2, c='coral', alpha=0.5, label=label2, edgecolor='coral')
    ax1.set_xlabel("Feature Score")
    ax1.set_ylabel("Feature reads/SP total")
    ax1.legend()
    ax1.set_title(title)
    plt.show()

def get_scores(merged_table):
    FiveP_scores = []
    ThreeP_scores = []
    j = 0
    for name in merged_table.columns:
        #print name
        if name.split(" ")[-1] =="Score":
            if j == 0:
                for score in merged_table[name]:
                    FiveP_scores.append(score)
                j += 1
            elif j == 1:
                for score in merged_table[name]:
                    ThreeP_scores.append(score)
                j += 1
            else:
                break
    both_scores = [FiveP_scores, ThreeP_scores]
    return both_scores

scores_df = SPTools.build_intron_tables(sys.argv[1])
exon_merged_df = merge_tables(2)
intron_merged_df = merge_tables(3)

exon_FiveP_scores = get_scores(exon_merged_df)[0]
exon_ThreeP_scores = get_scores(exon_merged_df)[1]
exons1 = get_points(exon_merged_df, 0)
exons2 = get_points(exon_merged_df, 1)
intron_FiveP_scores = get_scores(intron_merged_df)[0]
intron_ThreeP_scores = get_scores(intron_merged_df)[1]
introns1 = get_points(intron_merged_df, 0)
introns2 = get_points(intron_merged_df, 1)

#print len(exon_FiveP_scores)
#print len(exons1)

lengths_ex= []
for name in exon_merged_df.columns:
    if name == "Intron Length":
        for entry in exon_merged_df[name]:
               lengths_ex.append(entry)
#print len(lengths_ex)

lengths_int=[]
for name in intron_merged_df.columns:
    if name == "Intron Length":
        for entry in intron_merged_df[name]:
               lengths_int.append(entry)
#print len(lengths_int)

scatter_plot(exon_FiveP_scores, exons1, exons2, "Replicate 1", "Replicate 2", "5'ss scores vs. normalized 5' exon reads")
scatter_plot(intron_FiveP_scores, introns1, introns2, "Replicate 2", "Replicate 2", "5'ss scores vs. normalized intron reads")
scatter_plot(lengths_ex, exons1, exons2, "Replicate 1", "Replicate 2", "Intron lengths vs. normalized 5' exon reads")
scatter_plot(exon_ThreeP_scores, exons1, exons2, "Replicate 1", "Replicate 2", "3'ss scores vs. normalized 5' exon reads")
scatter_plot(intron_ThreeP_scores, introns1, introns2, "Replicate 1", "Replicate 2", "3'ss scores vs. normalized intron reads")
scatter_plot(lengths_int, introns1, introns2, "Replicate 1", "Replicate 2", "Intron lengths vs. normalized intron reads")