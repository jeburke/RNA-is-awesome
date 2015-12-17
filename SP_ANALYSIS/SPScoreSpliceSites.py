__author__ = 'jordanburke'

''' This script compares splice site scores and normalized tables from SPNormalizeAndPlot and plots the results
Usage: NormalizeAndPlot.py <splice site scores> <normalized exon table> <normalized intron table> <prefix> '''

import sys
import SPTools
import pandas as pd
import matplotlib.pyplot as plt

n = 2
scores_df = SPTools.build_intron_tables(sys.argv[1])

while n < len(sys.argv)-1:
    file = sys.argv[n]
    counts_df = SPTools.build_intron_tables(file)
    merged_df = pd.merge(scores_df, counts_df, right_index=True, left_index=True)
    n += 1
scores_df_sorted = merged_df.sort(columns=scores_df.columns[0])
fout = open("{0}_scored_ss.txt".format(sys.argv[-1].split(".")[0]), "w")
print fout
scores_df_sorted.to_csv(fout, sep='\t')
fout.close()

def get_points(normalizedResults, a):
    values = []
    n = 0
    for name in normalizedResults.columns:
        if name.split("_")[-1] == "norm":
            print name
            if n == a:
                s = pd.Series(normalizedResults[name])
                values = s.tolist()
                n += 1
            else:
                n += 1
    return values

FiveP_scores = []
ThreeP_scores = []
j = 0
for name in scores_df_sorted.columns:
    if name.split("_")[-1] =="Score":
        #print scores_df_sorted[name]
        #print name
        if j == 0:
            #s = pd.Series(scores_df_sorted[name])
            #FiveP_scores = s.tolist()
            for score in scores_df_sorted[name].iteritems():
                FiveP_scores.append(score)
            #print FiveP_scores
            j += 1
        elif j == 1:
            #s = pd.Series(scores_df_sorted[name])
            #print s
            for score in scores_df_sorted[name].iteritems():
                ThreeP_scores.append(score)
            #ThreeP_scores = s.toList()
            j += 1
        else:
            break

exons1 = get_points(scores_df_sorted, 0)
exons2 = get_points(scores_df_sorted, 1)

print FiveP_scores[0:20]
print len(ThreeP_scores)
print exons1[0:20]
print len(exons2)


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.scatter(FiveP_scores, exons1, c='royalblue', label="5' splice site scores", alpha = 0.5, edgecolor='darkslateblue')
ax1.scatter(ThreeP_scores, exons1, c='coral', alpha=0.5, label="3' splice site scores", edgecolor='coral')
ax1.set_xlabel("Splice site scores")
ax1.set_ylabel("Replicate 2 (log10)")
#xlim = ax1.get_xlim()
#ax1.set_ylim(xlim)
#ax1.set_xlim(xlim)
#ylim = ax1.get_ylim()
ax1.legend()
ax1.set_title("Splice site scores vs 5' exon values")
#ax1.plot(xlim, ylim, ls="--", c=".3", lw=1.5)
plt.show()