import sys

def build_transcript_dic(gff):
    roman_num = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"]
    latin_num = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16"]
    num_dict = dict(zip(roman_num, latin_num))

    feature_list = ["CDS","intron","exon"]
    counter = 0
    isoform_di={}
    with open(gff,"r") as fin:
        for line in fin:
            if line.startswith('chr'):
                columns = line.split("\t")
                line_list = []
                if columns[0][3:].strip() in roman_num:
                    number = num_dict[columns[0][3:]]
                    new_chromosome = "chr"+number
                elif columns[0][3:].strip() in latin_num:
                    new_chromosome = columns[0].strip()
                
                gene_ID = columns[8].strip()
                if columns[2].strip() == "mRNA":
                    gene_ID = gene_ID.split(";")[0]
                    gene_ID = gene_ID.split("=")[1]
                    if gene_ID.endswith("mRNA"):
                        gene_ID = gene_ID.split("_")[0]
                    if gene_ID not in isoform_di:
                        isoform_di[gene_ID] = []
                    line_list.append(new_chromosome)
                    line_list = line_list+columns[1:9]
                    isoform_di[gene_ID].append(line_list)
                          
                elif columns[2].strip() in feature_list:
                    if columns[2].strip() == "intron":
                        counter += 1
                    if gene_ID.startswith("Parent"):
                        gene_ID = gene_ID.split(";")[0]
                        gene_ID = gene_ID.split("=")[1]
                    elif gene_ID.startswith("ID"):
                        gene_ID = gene_ID.split(";")[1]
                        gene_ID = gene_ID.split("=")[1]
 
                    line_list.append(new_chromosome)
                    line_list = line_list+columns[1:9]
                    if "CNAG" not in gene_ID:
                        if gene_ID.split("_")[0] not in isoform_di:
                            isoform_di[gene_ID.split("_")[0]] = []
                        isoform_di[gene_ID.split("_")[0]].append(line_list)
                    elif "CNAG" in gene_ID:
                        if gene_ID not in isoform_di:
                            isoform_di[gene_ID] = []
                        isoform_di[gene_ID].append(line_list)

    print counter
    return isoform_di

def get_transcript_lengths(isoform_di, prefix):
    length_dict = {}
    for gene_ID, lines in isoform_di.iteritems():
        for line in lines:
            if line[2] == "mRNA":
                length_dict[gene_ID] = str(int(line[4])-int(line[3]))
    with open("{0}_transcript_lengths.txt".format(prefix), "w") as fout:
        fout.write("Transcript\t Length\n")
        for gene_ID, length in length_dict.iteritems():
            fout.write(gene_ID+"\t"+length+"\n")
            
def get_intron_lengths(isoform_di, prefix):
    length_dict = {}
    counter = 0
    for gene_ID, lines in isoform_di.iteritems():
        exon_list = []
        intron_list = []
        for line in lines:
            if line[2] == "intron":
                counter += 1
                intron_list.append(str(int(line[4])-int(line[3])))
                
            elif line[2] == "exon":
                exon_list.append(int(line[3]))
                exon_list.append(int(line[4]))
        if exon_list > 2:
            exon_list.sort()
            n = 2
            while n < len(exon_list):
                intron_list.append(exon_list[n]-exon_list[n-1])
                n += 2
        if len(intron_list) > 0:
            if line[6] == "-":
                intron_list.reverse()
            length_dict[gene_ID] = intron_list
        
    print counter
    with open("{0}_intron_lengths.txt".format(prefix), "w") as fout:
        fout.write("Transcript\t Intron\t Length\n")
        for gene_ID, lengths in length_dict.iteritems():
            n = 0
            while n < len(lengths):
                line_list = [gene_ID, str(n+1), str(lengths[n]),"\n"]
                fout.write("\t".join(line_list))
                n += 1
            
                
isoform_di = build_transcript_dic(sys.argv[1])
get_transcript_lengths(isoform_di, sys.argv[2])
get_intron_lengths(isoform_di, sys.argv[2])