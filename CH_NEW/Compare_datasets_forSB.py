'''Script to graph various metrics for comparing datasets
python Compare_datasets.py telomere_summaries_for_datasets telomere_summary_for_WCE prefix_for_wildtype

'''

import sys
import matplotlib.pyplot as plt
import numpy as np
import math

centromeres = {}
centromeres["I"] = [3753687, 3789421]
centromeres["II"] = [1602264, 1644747]
centromeres["III"] = [1070904, 1137003]

total_cen_len = 0
for centromere in centromeres:
    cen_data = centromeres[centromere]
    length = cen_data[1] - cen_data[0]
    total_cen_len += length

#Read in WCE data
wce_tel_5P_bg = {}
wce_tel_3P_bg = {}
wce_cen_bg = {}
genome = {}
wce_ave_coverage = 0
wce_regions = 0

#Read in Tel Summary file first
fin = open(sys.argv[-1], "r")
fin.readline()
for line in fin:
    data = line.split("\t")
    if len(data) > 7:
        if data[1] == "3'":
            genome[data[0]] = int(data[3])
        ###Not sure what this line is doing
        #wce_ave_coverage = float(data[5])
        wce_ave_coverage += float(data[5])
        wce_regions += 1
    elif len(data) == 4:
        wce_telo_ness = float(data[0].split(" ")[1])
        wce_cen_ness = float(data[1].split(" ")[1])
        wce_total_reads = int(data[2].split(" ")[2])
            
fin.close()

wce_ave_coverage = float(wce_ave_coverage)/wce_regions


#Now read in coverage data from WCE
for chrom_name in centromeres:
    wce_tel_5P = np.zeros([1,80000])[0]
    wce_tel_3P = np.zeros([1,80000])[0]
    len_cen = int(centromeres[chrom_name][1] - centromeres[chrom_name][0])
    wce_cen = np.zeros([1,len_cen])[0]
    
    wce_tel_5P_bg[chrom_name] = wce_tel_5P
    wce_tel_3P_bg[chrom_name] = wce_tel_3P
    wce_cen_bg[chrom_name] = wce_cen
    
    
print "Reading in WCE1"
wce1_fin = open("{0}.bedgraph".format(sys.argv[-1].split("_telomere_")[0]), "r")
for line in wce1_fin:
    data = line.split("\t")
    chrom_name = data[0]
    if chrom_name == "MT" or chrom_name == "MTR" or chrom_name == "AB325691":
        pass
    else:
        start_3prime = int(genome[chrom_name]) - 80000
        start = int(data[1])
        end = int(data[2])
        coverage = float(data[3])
        centromere = centromeres[chrom_name]
        cen_center = int(centromere[0] + ((centromere[1] - centromere[0]) / 2.))
        if end < 80000:
            wce_tel_5P_bg[chrom_name][start] += coverage
        elif end > start_3prime:
            wce_tel_3P_bg[chrom_name][int(genome[chrom_name]) - end ] += coverage
        elif end > centromere[0] and end <= centromere[1]:
            wce_cen_bg[chrom_name][end - centromere[0] - 1 ] += coverage
wce1_fin.close()

telo_list = []
telo_list_rpkm = []
cen_list = []
cen_list_rpkm = []
labels = []

data_array_5P = {}
data_array_3P = {}
data_array_5P_scaled = {}
data_array_3P_scaled = {}
data_array_cen_scaled = {}
data_array_cen = {}

for m, arg in enumerate(sys.argv[1:]):
    print m
    tel_bgs_5P = {}
    tel_bgs_3P = {}
    tel_bgs_5P_scaled = {}
    tel_bgs_3P_scaled = {}
    cen_bg = {}
    cen_bg_scaled = {}
    for chrom in centromeres:
        tel_bgs_5P[chrom] = np.zeros([1,80000])[0]
        tel_bgs_3P[chrom] = np.zeros([1,80000])[0]
        tel_bgs_5P_scaled[chrom] = np.zeros([1,80000])[0]
        tel_bgs_3P_scaled[chrom] = np.zeros([1,80000])[0]
        
        len_cen = centromeres[chrom][1] - centromeres[chrom][0]
        cen_bg[chrom] = np.zeros([1,len_cen])[0]
        cen_bg_scaled[chrom] = np.zeros([1,len_cen])[0]

    #Get average coverage
    fin = open(arg, "r")
    genome = {}

    name_string = arg.split("/")[-1].split("_")
    genotype = name_string[0]
    for name in name_string[1:]:
        if name != "perbase":
            genotype = genotype + "_" + name
        else:
            break
    print genotype
    
    if m  == 0:
        wildtype_genotype = genotype
        print "Wildtype: {0}".format(wildtype_genotype)
    
    if m == len(sys.argv[1:])-1:
        wce_genotype = genotype
        print "WCE: {0}".format(wce_genotype)
    
    fin.readline()
    total_tel_len = 0
    ave_coverage = 0
    line_count = 0
    for line in fin:
        data = line.split("\t")
        if len(data) > 7:
            if data[1] == "3'":
                genome[data[0]] = int(data[3])
            #ave_coverage = float(data[5])
            ave_coverage += float(data[5])
            line_count += 1
            if data[4] != "":
                total_tel_len += float(data[4])
        elif len(data) == 4:
            telo_ness = float(data[0].split(" ")[1])
            cen_ness = float(data[1].split(" ")[1])
            total_reads = float(data[2].split(" ")[2])
    ave_coverage = float(ave_coverage)/line_count
            
    if telo_ness > 0:
        telo_list.append(telo_ness / total_reads * 100000)
    else:
        telo_list.append(0)
        
    if cen_ness > 0:
        cen_list.append(cen_ness / total_reads * 100000)
    else:
        cen_list.append(0)

    labels.append(genotype)
    fin.close()
    
    fin = open("{0}.bedgraph".format("_".join(arg.split("_")[:-2])))
    norm = wce_ave_coverage / ave_coverage
    for line in fin:
        data = line.split("\t")
        chrom_name = data[0]
        if chrom_name == "MT" or chrom_name == "MTR" or chrom_name == "AB325691":
            pass
        else:
            start_3prime = int(genome[chrom_name]) - 80000
            start = int(data[1])
            end = int(data[2])
            coverage = float(data[3])
            centromere = centromeres[chrom_name]
            cen_center = int(centromere[0]+((centromere[1] - centromere[0]) / 2.))
            if end < 80000:
                tel_bgs_5P[chrom_name][start] += coverage
                tel_bgs_5P_scaled[chrom_name][start] += (coverage / wce_tel_5P_bg[chrom_name][start] * norm)
            elif end > start_3prime:
                tel_bgs_3P[chrom_name][int(genome[chrom_name]) - end ] += coverage
                tel_bgs_3P_scaled[chrom_name][int(genome[chrom_name]) - end] += (coverage / wce_tel_3P_bg[chrom_name][int(genome[chrom_name]) - end] * norm)
            elif end > centromere[0] and end <= centromere[1]:
                cen_bg[chrom_name][end - centromere[0] - 1 ] += coverage
                cen_bg_scaled[chrom_name][end - centromere[0] - 1] += (coverage / wce_cen_bg[chrom_name][end - centromere[0] - 1] * norm)

    data_array_5P_scaled[genotype] = tel_bgs_5P_scaled
    data_array_3P_scaled[genotype] = tel_bgs_3P_scaled
    data_array_cen_scaled[genotype] = cen_bg_scaled
    fin.close()

colors = {}
zorders = {}
alphas = {}
mut_colors = ['red','orange','gold', 'yellow']
print "Number of samples: "+str(len(data_array_3P_scaled))
n=0
print data_array_3P_scaled.keys()
for key in data_array_3P_scaled:
    if key == wce_genotype:
        pass
    elif key == wildtype_genotype:
        print "wildtype is detected"
        colors[key] = '0.2'
        alphas[key] = 1.0
        zorders[key] = 1
    else:
        colors[key] = mut_colors[n]
        zorders[key] = 2
        alphas[key] = 1
        n += 1
        
#print len(data_array_3P_scaled)
#print data_array_3P_scaled

total_tel = {}
total_cen = {}
for chrom in centromeres:
    tel_wce_5P_bg_scaled = []
    tel_wce_3P_bg_scaled = []
    cen_wce_bg_scaled = []

    save_name_list = []
    all_y = []
    for genotype in data_array_5P_scaled:
        if genotype != wce_genotype:
            print genotype
            save_name_list.append(genotype)
            tel_bedgraph_scaled = data_array_5P_scaled[genotype][chrom]
            ind = np.arange(len(tel_bedgraph_scaled))
            ind = [x/1000. for x in ind]
            plt.plot(ind, tel_bedgraph_scaled, color=colors[genotype], zorder= zorders[genotype], alpha=alphas[genotype])
            all_y.append(max(tel_bedgraph_scaled))
            if chrom != 'III':
                if genotype not in total_tel:
                    total_tel[genotype] = 0
                total_tel[genotype] += sum(tel_bedgraph_scaled)/float(len(tel_bedgraph_scaled))
        else:
            tel_wce_5P_bg_scaled = data_array_5P_scaled[genotype][chrom]
    save_name_prefix = "_vs_".join(save_name_list)
    plt.xlim(0,len(tel_bedgraph_scaled)/1000.)
    plt.ylim(0.5, max(all_y)+0.1*max(all_y))
    #plt.ylim(1,max(tel_bedgraph_scaled))
    save_name = "Telomere_5P_chrom{0}_coverage_scaled_{1}.pdf".format(chrom,save_name_prefix)
    plt.savefig(save_name, close=True, verbose=True, format='pdf')
    plt.title(chrom+' Left Telomere')
    plt.xlabel('Distance from telomere (kb)')
    plt.ylabel('ChIP signal')
    plt.show()
    plt.clf()
    
    #ind = np.arange(len(tel_wce_5P_bg_scaled ))
    #plt.bar(ind, tel_wce_5P_bg_scaled, width = 1.0, color='grey', linewidth = 0)
    #plt.xlim(0,len(tel_wce_5P_bg_scaled))
    #save_name = "Telomere_5P_chrom{0}_wce_coverage_scaled_{1}".format(chrom,save_name_prefix)
    #plt.savefig(save_name, close=True, verbose=True)
    #plt.clf()

    all_y = []
    for genotype in data_array_3P_scaled:
        if genotype != wce_genotype:
            print genotype
            tel_bedgraph_scaled = data_array_3P_scaled[genotype][chrom]
            ind = np.arange(len(tel_bedgraph_scaled))
            ind = [x/1000. for x in ind]
            plt.plot(ind, tel_bedgraph_scaled, color=colors[genotype], zorder= zorders[genotype], alpha=alphas[genotype])
            all_y.append(max(tel_bedgraph_scaled))
            if chrom != 'III':
                if genotype not in total_tel:
                    total_tel[genotype] = 0
                total_tel[genotype] += sum(tel_bedgraph_scaled)/float(len(tel_bedgraph_scaled))
        else:
            tel_wce_5P_bg_scaled = data_array_5P_scaled[genotype][chrom]
    save_name = "Telomere_3P_chrom{0}_coverage_scaled_{1}.pdf".format(chrom,save_name_prefix)
    plt.xlim(0,len(tel_bedgraph_scaled)/1000.)
    plt.ylim(0.5, max(all_y)+0.1*max(all_y))
    plt.savefig(save_name, close=True, verbose=True, format='pdf')
    plt.title(chrom+' Right Telomere')
    plt.xlabel('Distance from telomere (kb)')
    plt.ylabel('ChIP signal')
    plt.show()
    plt.clf()

    #ind = np.arange(len(tel_wce_3P_bg_scaled ))
    #plt.bar(ind, tel_wce_3P_bg_scaled, width = 1.0, color='grey', linewidth = 0)
    #save_name = "Telomere_3P_chrom{0}_wce_coverage_scaled_{1}".format(chrom,save_name_prefix)
    #plt.savefig(save_name, close=True, verbose=True)
    #plt.clf()

    all_y = []
    for genotype in data_array_cen_scaled:
        if genotype != wce_genotype:
            print genotype
            cen_bedgraph_scaled = data_array_cen_scaled[genotype][chrom]
            #ind = np.arange(len(cen_bedgraph_scaled))
            ind = range(0-len(cen_bedgraph_scaled)/2,len(cen_bedgraph_scaled)/2)
            if len(ind) != len(cen_bedgraph_scaled):
                ind.append(max(ind)+1)
            ind = [x/1000. for x in ind]
            plt.plot(ind, cen_bedgraph_scaled, color=colors[genotype], zorder= zorders[genotype], alpha=alphas[genotype])
            all_y.append(max(cen_bedgraph_scaled))
            if genotype not in total_cen:
                total_cen[genotype] = 0
            total_cen[genotype] += sum(cen_bedgraph_scaled)/float(len(cen_bedgraph_scaled))
        else:
            cen_wce_5P_bg_scaled = data_array_cen_scaled[genotype][chrom]
    save_name = "Centromere_chrom{0}_coverage_scaled_{1}.pdf".format(chrom, save_name_prefix)
    plt.xlim(min(ind),max(ind))
    plt.savefig(save_name, close=True, verbose=True, format='pdf')
    plt.ylim(0.5, max(all_y)+0.1*max(all_y))
    plt.title(chrom+' Centromere')
    plt.xlabel('Distance from centromere (kb)')
    plt.ylabel('ChIP signal')
    plt.show()
    plt.clf()
 
    #ind = np.arange(len(cen_wce_bg_scaled ))
    #plt.bar(ind, cen_wce_bg_scaled, width = 1.0, color='grey', linewidth = 0)
    #plt.xlim(0,len(cen_wce_bg_scaled))
    #save_name = "Centromere_chrom{0}_wce_coverage_scaled_{1}".format(chrom,save_name_prefix)
    #plt.savefig(save_name, close=True, verbose=True)
    #plt.clf()
    
#for term in telo_list:
#    telo_list_rpkm.append(term / total_tel_len)
#for term in cen_list:
#    cen_list_rpkm.append(term / total_cen_len)



#print "Telo list: {0}".format(telo_list)
#print "Cen list: {0}".format(cen_list)
width = 0.5

#telo_list = telo_list[:-1]
#cen_list = cen_list[:-1]
#telo_list_rpkm = telo_list_rpkm[:-1]
#cen_list_rpkm  = cen_list_rpkm[:-1]
#labels=labels[:-1]
labels = []
telo_list = []
cen_list = []
for key, value in total_tel.iteritems():
    labels.append(key)
    telo_list.append(value)
    cen_list.append(total_cen[key])

ind = np.arange(len(telo_list))
fig, axes = plt.subplots(nrows=2, sharex=True)
y_max = 1.1*max(telo_list+cen_list)
axes[0].bar(ind, telo_list, align='center', color='navy')
axes[1].bar(ind, cen_list, align='center', color='gold')
axes[0].set(xticks=range(len(labels)), xticklabels=labels)
axes[1].set(xticks=range(len(labels)), xticklabels=labels)
#axes[0].xaxis.tick_bottom()
#axes[1].xaxis.tick_top()
axes[0].set_ylim(0,y_max)
axes[1].set_ylim(0,y_max)
fig.subplots_adjust(hspace=0)
axes[1].invert_yaxis()
save_name = "Compare_all_{0}.svg".format(save_name_prefix)
plt.savefig(save_name, format="svg", close=True, verbose = True)
plt.show()
plt.clf()


#ind = np.arange(len(telo_list_rpkm))
#fig, axes = plt.subplots(nrows=2, sharex=True)
#axes[0].bar(ind, telo_list_rpkm, align='center', color='navy')
#axes[1].bar(ind, cen_list_rpkm, align='center', color='gold')
#axes[1].invert_yaxis()
#axes[0].set(xticks=range(len(labels)), xticklabels=labels)
#axes[0].xaxis.tick_bottom()
#axes[1].xaxis.tick_top()
#fig.subplots_adjust(hspace=0.3)
#save_name = "Compare_all_rpkm_{0}.svg".format(save_name_prefix)
#plt.savefig(save_name, format="svg", close=True, verbose = True)
#plt.show()
#plt.clf()