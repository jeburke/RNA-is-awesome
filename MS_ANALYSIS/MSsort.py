__author__ = 'jordanburke'
import csv
import sys


f_data = sys.argv[1]
f_list = sys.argv[2]
fout = ("{0}_sorted.tsv".format(sys.argv[1].split(".")[0]))


data_dict = {}
with open(f_data, 'rU') as data:
    datareader = csv.reader(data, delimiter="\t")
    counter = 1
    for row in datareader:
        if counter == 1:
            with open(fout, 'w') as f:
                f.write("Name\t"+"\t".join(row)+"\n")
            counter +=1

        data_dict[row[0].strip()] = []
        n = 1
        while n < len(row):
            data_dict[row[0].strip()].append(row[n])
            n +=1

new_list = []
with open(f_list, 'rU') as gene_list:
    listreader = csv.reader(gene_list, delimiter="\t")
    for row in listreader:
        new_list.append(row[0].strip())


for gene in new_list:
    with open(fout, 'a') as f:
        line_list = []
        line_list.append(gene)
        line_list.append(gene)
        for i in range(len(data_dict[gene])):
            line_list.append(data_dict[gene][i])
            i += 1
        line_list.append("\n")
        line = "\t".join(line_list)
        f.write(line)



