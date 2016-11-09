import sys
import os
import csv

def reformat_file(sf_file, suffix_list):
    new_line_dict = {}
    with open(sf_file, 'r') as fin:
	n=0
	for line in fin:
	    if n == 0: n += 1
	    else:
		columns = line.split('\t')
		name = columns[0]
		num_reads = int(float(columns[4].strip()))
		if len(suffix_list) > 0:
		    flag = False
		    for suffix in suffix_list:
			if suffix not in new_line_dict:
			    new_line_dict[suffix] = []
			
			if suffix in name:
			    name = name.split('_'+suffix)[0]
			    new_line_dict[suffix].append([name,str(num_reads),'\n'])
			    flag = True
		    if flag == False:
			if 'all' not in new_line_dict: new_line_dict['all'] = []
			new_line_dict['all'].append([name,str(num_reads),'\n'])
		else:
		    if 'all' not in new_line_dict: new_line_dict['all'] = []
		    new_line_dict['all'].append([name,str(num_reads),'\n']) 
    return new_line_dict

def write_new_files(sf_file, new_line_dict):
    for suffix, lines in new_line_dict.iteritems():
	with open('{0}{1}_{2}.tsv'.format(directory, sf_file.split('/')[-1].split('.')[0], suffix), 'w') as fout:
	    for line in lines:
		fout.write('\t'.join(line))

directory = sys.argv[1]
file_list = []

for file in os.listdir(directory):
    if file.endswith('sf'):
        file_list.append(directory+file)

suffix_list = []
if len(sys.argv) > 2:
   with open(sys.argv[2], 'r') as suffixes:
	for line in suffixes:
	    suffix_list.append(line.strip())

for sf_file in file_list:
    print sf_file
    new_line_dict = reformat_file(sf_file, suffix_list)    
    write_new_files(sf_file, new_line_dict)
