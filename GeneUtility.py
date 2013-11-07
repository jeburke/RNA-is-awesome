'''
Created on Aug 19, 2012
Written using Python 2.7.3
@author: chomer
'''
CNAG_index = 0
Gene_start_index = 4
Gene_end_index = 5
Gene_strand_index = 6
Gene_chromosome_index = 8
Peak_length_index = 3
Abs_summit_index = 4
Pileup_index = 5
Peak_chromosome_index = 0
Log_p_value_index = 6
Fold_enrichment_index = 7
Log_q_value_index = 8
peak_name_index = 9
CNAG_annotation_index = 0
gene_name_index = 7
pfam_numbers_index = 16
chromend=[2291499, 1621675, 1575141, 1084805, 1814975, 1422463, 1399503, 1398693, 1186808, 1059964, 1561994, 774062, 756744, 926563]    
import sys
from pylab import *


class CommonEquality(object):
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False
    
    def __ne__(self,other):
        return not self.__eq__(other)
    
class GeneAnnotationDict(CommonEquality):
    '''Dictionary of Gene Annotations so we can look up annotations by CNAG'''
    def __init__(self, AnnotationsByCNAG={}):
        self.AnnotationsByCNAG = AnnotationsByCNAG

    def __iter__(self):
        return iter(self.AnnotationsByCNAG)
    
    def PopulateGeneAnnotation_broad(self):
        filename = "/home/chomer/Code_Base/PeakLocator/cryptococcus_neoformans_grubii_h99_2_genome_summary_per_gene.txt"
        with open(filename, "r") as f:
            f.next()
            for line in f:
                data = line.split("\t")
                CNAG = data[CNAG_annotation_index]
                gene_name = data[gene_name_index]
                pfam_numbers = data[pfam_numbers_index]
                self.AnnotationsByCNAG[CNAG] = GeneAnnotation(CNAG, gene_name, pfam_numbers)
        f.close()
                
    def PopulateGeneAnnotation_new(self):
        filename1 = "/home/chomer/Code_Base/PrimerDesigner/CNA2_FINAL_CALLGENES_2.gff3"
        filename2 = "/home/chomer/Code_Base/PrimerDesigner/CNA2_MISC_RNA_1.gff3"
        with open(filename1, "r") as f:
            for line in f:
                data = line.split("\t")
                if len(data) > 1:
                    if data[2] == "gene":
                        CNAG = data[8].split(";")[0].split("=")[1]
                        gene_name = data[8].split(";")[2].split("=")[1].rstrip()
                        pfam_numbers = ""
                        self.AnnotationsByCNAG[CNAG] = GeneAnnotation(CNAG, gene_name, pfam_numbers)
        f.close()
        with open(filename2, "r") as f:
            for line in f:
                data = line.split("\t")
                if len(data) > 1:
                    if data[2] == "gene":
                        CNAG = data[8].split(";")[0].split("=")[1]
                        gene_name = data[8].split(";")[2].split("=")[1].rstrip()
                        pfam_numbers = ""
                        self.AnnotationsByCNAG[CNAG] = GeneAnnotation(CNAG, gene_name, pfam_numbers)
                
    def LookUpGeneAnnotationByCNAG(self, CNAG):
        return self.AnnotationsByCNAG[CNAG]     

class GeneAnnotation(CommonEquality):
    '''Gene Annotation object that makes up above dictionary'''
    def __init__(self, CNAG, gene_name, pfam_numbers):
        self.CNAG = CNAG
        self.gene_name = gene_name
        self.pfam_numbers = pfam_numbers
        
    def __str__(self):
        '''for printing object as string'''
        return '%s %s %s' % (self.CNAG, self.gene_name, self.pfam_numbers)

class IntergenicDict(CommonEquality):
    '''Dictionary of Intergenic Objects and functions for manipulating them'''
    def __init__(self, IntergenicsByCNAG={}):
        self.IntergenicsByCNAG = IntergenicsByCNAG
        
    def __iter__(self):
        return iter(self.IntergenicsByCNAG)
        
    def FindIntergenicByCNAG(self,CNAG):
        return self.IntergenicsByCNAG[CNAG]

    def PopulateIntergenicRegion(self, genedict):
        genes = genedict.GetGenes()
        CNAGs = []
        intergenic_sizes = []
        for gene in genes: 
            CNAGs.append(gene.CNAG)
        for CNAG in CNAGs:
            strand = genedict.FindGeneStrand(CNAG)
            intergenic = genedict.FindIntergenicRegion(CNAG, strand)
            self.IntergenicsByCNAG[CNAG]=intergenic
            length = int(intergenic.start) - int(intergenic.end)
            intergenic_sizes.append(abs(length))
        print "Average size of intergenics in crypto = {0}".format(sum(intergenic_sizes)/float(len(intergenic_sizes)))
        
    def PopulateNonPromoterIntergenicRegions(self, genedict):
        genes = genedict.GetGenes()
        for gene in genes: 
            strand = gene.strand
            neighbor_gene = genedict.FindNeighboringGene(gene.CNAG, "downstream")
#            if neighbor_gene != None:
#                print "gene: {0} downstream neighbor: {1}".format(gene.CNAG, neighbor_gene.CNAG)
            if neighbor_gene != None:
                if gene.strand == "+" and gene.end < neighbor_gene.start and neighbor_gene.strand == "-":
                    intergenic1 = Intergenic("3P{0}".format(neighbor_gene.CNAG), neighbor_gene.start, gene.end, neighbor_gene.strand, neighbor_gene.chromosome)
                    intergenic2 = Intergenic("3P{0}".format(gene.CNAG), gene.end, neighbor_gene.start, gene.strand, gene.chromosome  )
                    self.IntergenicsByCNAG["3P{0}".format(neighbor_gene.CNAG)] = intergenic1
                    self.IntergenicsByCNAG["3P{0}".format(gene.CNAG)] = intergenic2
                elif gene.strand == "-" and gene.start > neighbor_gene.end and neighbor_gene.strand == "+":
                    intergenic1 = Intergenic("3P{0}".format(neighbor_gene.CNAG), neighbor_gene.end, gene.start, neighbor_gene.strand, neighbor_gene.chromosome)
                    intergenic2 = Intergenic("3P{0}".format(gene.CNAG), gene.start, neighbor_gene.end, gene.strand, gene.chromosome)
                    self.IntergenicsByCNAG["3P{0}".format(neighbor_gene.CNAG)] = intergenic1
                    self.IntergenicsByCNAG["3P{0}".format(gene.CNAG)] = intergenic2

    
    def Update(self, other_dict):
        self.IntergenicsByCNAG.update(other_dict)
        
    def FindTypes(self, gene_dict, misc_dict):
        '''Populate the type of each intergenic region depending on gene orientation and miscRNA orientation/position'''
        num_1 = 0
        num_2 = 0
        num_3 = 0
        list_of_types = []
        typecountdict = {}
        for CNAG in self.IntergenicsByCNAG:
            intergenic = self.FindIntergenicByCNAG(CNAG)
            '''First find the overall promoter type'''
            if intergenic.CNAG != None and intergenic.CNAG[0:4] == "CNAG" :
                upstream_gene = gene_dict.FindNeighboringGene(intergenic.CNAG, "upstream")
                if upstream_gene != None and upstream_gene.strand == intergenic.strand:
                    intergenic_type = "Type1"
                    num_1 += 1
                elif upstream_gene != None and upstream_gene.strand != intergenic.strand:
                    intergenic_type = "Type2"
                    num_2 += 1 
            elif intergenic.CNAG != None and intergenic.CNAG[0:4] == "3PCN":
                intergenic_type = "Type3"
                num_3 += 1
            '''Next find the Promoter Misc Type'''
            for miscRNA in misc_dict.FilterByChromosome(intergenic.chromosome):
                if intergenic.start > intergenic.end:
                    intergenic_temp = intergenic.start
                    intergenic.start = intergenic.end
                    intergenic.end = intergenic_temp
                if miscRNA.strand == "-":
                    misc_temp = miscRNA.start
                    miscRNA.start = miscRNA.end
                    miscRNA.end = misc_temp
                if (int(miscRNA.start) > int(intergenic.start)) and (int(miscRNA.end) < int(intergenic.end)) or (int(miscRNA.start) > (int(intergenic.start) - (abs(int(intergenic.start) - int(intergenic.end))*0.20)) and (int(miscRNA.end) < int(intergenic.end))) or ((int(miscRNA.start) > int(intergenic.start)) and int(miscRNA.end) < (int(intergenic.end) + (abs(int(intergenic.start) - int(intergenic.end))*0.20))):
                    intergenic_type = intergenic_type+"a"
                elif ((int(miscRNA.start) > int(intergenic.start) - 850 and int(miscRNA.end) < int(intergenic.start) + (abs(int(intergenic.start) - int(intergenic.end))*0.20)) and miscRNA.strand == "-"):
                    intergenic_type = intergenic_type+"b"
                elif ((int(miscRNA.end) < int(intergenic.end) + 850 and int(miscRNA.start) > int(intergenic.end) - (abs(int(intergenic.start) - int(intergenic.end))*0.20)) and miscRNA.strand == "+"):
                    intergenic_type = intergenic_type+"b"
                elif ((int(miscRNA.start) > int(intergenic.start) - 850 and int(miscRNA.end) < int(intergenic.start) + (abs(int(intergenic.start) - int(intergenic.end))*0.20)) and miscRNA.strand == "+"):
                    intergenic_type = intergenic_type+"c"
                elif ((int(miscRNA.end) < int(intergenic.end) + 850 and int(miscRNA.start) > int(intergenic.end) - (abs(int(intergenic.start) - int(intergenic.end))*0.20)) and miscRNA.strand == "-"):
                    intergenic_type = intergenic_type+"c"
            list_of_types.append(intergenic_type)
            if intergenic_type in typecountdict:
                count = typecountdict[intergenic_type]
                if intergenic_type == "Type2" or intergenic_type == "Type3":
                    count = count + 0.5
                else:
                    count = count + 1
                typecountdict[intergenic_type] = count
            else:
                if intergenic_type == "Type2" or intergenic_type == "Type3":
                    typecountdict[intergenic_type] = 0.5
                else:
                    typecountdict[intergenic_type] = 1
            self.IntergenicsByCNAG[CNAG].type = intergenic_type
#        print "Number of Type 1 Intergenics: {0}".format(num_1)
#        print "Number of Type 2 Intergenics: {0}".format(num_2/2)
#        print "Number of Type 3 Intergenics: {0}".format(num_3/2)
#        print "Types of Promoters Found: {0}".format(typecountdict.keys())
#        print "Numbers of Each: {0}".format(typecountdict.values())
        
        '''Graph overall distribution of genome features in the genome'''
#        figure(1, figsize=(14,14))
#        ax = axes([0.1, 0.1, 0.8, 0.8])
 
#        labels = typecountdict.keys()
#        fracs = typecountdict.values()
#        def my_autopct(pct):
#            total=sum(fracs)
#            val=pct*total/100.0
#            return '{p:.0f}%  ({v:.0f})'.format(p=pct,v=val)
#        pie(fracs, labels=labels)
#        pie(fracs, labels=labels, autopct=my_autopct, shadow=True)
#        title('Promoter Types in C. neoformans Genome')
#        show()
#        savefig("Crypto_genome_promoter_distribution")
        '''Graph overall distribution of promoters in the genome'''
#        figure(11, figsize=(10,10)) 
#        ax = axes([0.1, 0.1, 0.8, 0.8])
#        labels = "Type1", "Type2", "Type3"
#        fracs = [num_1, num_2 / 2, num_3 / 2]
#        pie(fracs, labels=labels, autopct=my_autopct, shadow=True)
#        title('Promoter Types in Cryptococcus neoformans Genome')
#        savefig("Promoters_simple_Crypto")



class MiscDict(CommonEquality):
    '''List of misc RNA objects and functions for manipulating them'''
    def __init__(self, MiscsByCNAG={}):
        self.MiscsByCNAG = MiscsByCNAG

    def __getitem__(self, CNAG):
        return self.MiscsByCNAG[CNAG]
    
    def FilterByChromosome(self,chromosome):
        def MiscHasCorrectChromosome(miscRNA):
            return(chromosome==miscRNA.chromosome)
        misc_list=self.MiscsByCNAG.values()
        MiscWithCorrectChromosome=filter(MiscHasCorrectChromosome, misc_list)
        try: 
            return MiscWithCorrectChromosome
        except IndexError:
            print "There are no genes on chromosome %s" % (chromosome)
        
    def PopulateMiscFromFile(self):
        filename = "/home/chomer/Code_Base/PrimerDesigner/CNA2_MISC_RNA_1.gff3"
        start = sys.maxint
        end = 0
        CNAG = ""
        strand = ""
        chrom = ""
        with open(filename, "r") as f:
            for x, line in enumerate(f):
                data = line.split("\t")
                if len(data) > 1:
                    if data[2] == "gene":
                        if x != 0:
                            genetemp = Gene(CNAG, start, end, strand, chrom)
                            self.MiscsByCNAG[CNAG] = genetemp
                        CNAG = data[8].split(";")[0].split("=")[1]
                        strand = str(data[6])
                        chrom = data[0].split("_")[1][3:]
                        start = data[3]
                        end = data[4]
        genetemp = Gene(CNAG, start, end, strand, chrom)
        self.MiscsByCNAG[CNAG] = genetemp
        f.close()

    def PopulateIntergenics(self):
        misc_intergenic_dict = {}
        for miscRNA in self.GetMiscs():
            intergenic = self.FindIntergenicMisc(miscRNA)
            misc_intergenic_dict[miscRNA.CNAG] = intergenic
        return misc_intergenic_dict
    
    def __iter__(self):
        return iter(self.MiscsByCNAG)
    
    def GetDict(self):
        '''returns the dictionary'''
        return self.MiscsByCNAG
    
    def GetMiscs(self):
        '''returns the list of gene objects found in the dictionary'''
        #print "number = {0}".format(len(self.MiscsByCNAG.values()))
        return self.MiscsByCNAG.values()
        
    def FindIntergenicMisc(self, miscRNA):
        #A promoter region for misc RNAs is defined as 850 bp upstream of the start of the coding region because this is the average intergenic region size in crypto
        Chromosome = miscRNA.chromosome
        strand = miscRNA.strand
        if strand == "+":
            StartIntergenic = int(miscRNA.start) - 850
            EndIntergenic = miscRNA.start
        elif strand == "-":
            StartIntergenic=miscRNA.end
            EndIntergenic = int(miscRNA.end) + 850
        else:
            raise Exception("strand must be defined as '+' or '-' with the quotes")
        intergenic = Intergenic(miscRNA.CNAG, StartIntergenic, EndIntergenic, strand, Chromosome, "miscRNA promoter")
        return (intergenic)  

class GeneDict(CommonEquality):
    '''List of gene objects and functions for manipulating them'''
    def __init__(self, GenesByCNAG={}):
        self.GenesByCNAG = GenesByCNAG

    def __iter__(self):
        return iter(self.GenesByCNAG)
    
    def GetGenes(self):
        '''returns the list of gene objects found in the dictionary'''
        return self.GenesByCNAG.values()
    
    def PopulateFromFile_broad(self):
        '''designed to work with Broad Institute gene list file'''
        filename = "/home/chomer/Code_Base/PeakLocator/cryptococcus_neoformans_grubii_h99_2_genome_summary.txt"
        with open(filename, 'r') as f:
            f.next()                                            #to skip first line of file which is a header
            for line in f:
                data=line.split('\t')    
                genetemp = Gene(data[CNAG_index], data[Gene_start_index], data[Gene_end_index], data[Gene_strand_index], data[Gene_chromosome_index])
                self.GenesByCNAG[data[CNAG_index]]=genetemp
        f.close()
                
    def PopulateFromFile_new(self):
        '''designed to work with Guillem's new annotations'''
        filename = "/home/chomer/Code_Base/PrimerDesigner/CNA2_FINAL_CALLGENES_2.gff3"
        start = sys.maxint
        end = 0
        CNAG = ""
        strand = ""
        chrom = ""
        with open(filename, "r") as f:
            for x, line in enumerate(f):
                data = line.split("\t")
                if len(data) > 1:
                    if data[2] == "gene":
                        if x != 0:
                            genetemp = Gene(CNAG, start, end, strand, chrom)
                            self.GenesByCNAG[CNAG] = genetemp
                        start = sys.maxint
                        end = 0
                        CNAG = data[8].split(";")[0].split("=")[1]
                        strand = str(data[6])
                        chrom = data[0].split("_")[1][3:]
                    elif data[2] == "CDS":
                        ex_start = int(data[3])
                        ex_end = int(data[4])
                        if ex_start < start:
                            start = ex_start
                        if ex_end > end:
                            end = ex_end
        genetemp = Gene(CNAG, start, end, strand, chrom)
        self.GenesByCNAG[CNAG] = genetemp
        f.close()
        
        
    def __len__(self):
        return len(self.GenesByCNAG)

    def FindGeneLength(self,CNAG):
        gene=self.GenesByCNAG[CNAG]
        return len(gene)

    def FindGeneStrand(self, CNAG):
        gene=self.GenesByCNAG[CNAG]
        return gene.strand

    def FindGeneByCNAG(self,CNAG):
        return self.GenesByCNAG[CNAG]

    def FindGeneByStart(self,start,chromosome):
        def HasCorrectStartAndChromosome(gene):
            return(start==gene.start and chromosome==gene.chromosome)
        Gene_list=self.GenesByCNAG.values()
        GeneWithCorrectStartAndChromosome=filter(HasCorrectStartAndChromosome,Gene_list)
        try:
            return GeneWithCorrectStartAndChromosome[0]
        except IndexError:
            print "There is no gene starting at %s on chromosome %s" % (start, chromosome)
    
    def FilterByChromosome(self,chromosome):
        def HasCorrectChromosome(gene):
            return(chromosome==gene.chromosome)
        Gene_list=self.GenesByCNAG.values()
        GeneWithCorrectChromosome=filter(HasCorrectChromosome,Gene_list)
        try: 
            return GeneWithCorrectChromosome
        except IndexError:
            print "There are no genes on chromosome %s" % (chromosome)
            
    def FindNeighboringGene(self, CNAG, direction):
        '''This function takes a dictionary of Gene Objects, a particular CNAG, and the direction you want to find its nearest gene neighbor'''
        OriginalGene = self.FindGeneByCNAG(CNAG)
        Chromosome = OriginalGene.chromosome
        CandidateGenes = self.FilterByChromosome(Chromosome)
        NeighborDiff = sys.maxint
        NeighborLocation = ""
        if direction == "upstream":
            if OriginalGene.strand == "+":
                StartOri = OriginalGene.start
                for Gene in CandidateGenes:
                    StartNeighbor = Gene.start
                    Difference = StartOri - StartNeighbor
                    if Difference > 0 and Difference < NeighborDiff:
                        NeighborDiff = Difference
                        NeighborLocation = StartNeighbor
            elif OriginalGene.strand == "-":
                StartOri = OriginalGene.end
                for Gene in CandidateGenes:
                    StartNeighbor = Gene.start
                    Difference = StartNeighbor - StartOri
                    if Difference > 0 and Difference < NeighborDiff:
                        NeighborDiff = Difference
                        NeighborLocation = StartNeighbor
            
        elif direction == "downstream":
            if OriginalGene.strand == "+":
                StartOri=OriginalGene.end
                for Gene in CandidateGenes: 
                    StartNeighbor = Gene.start
                    Difference = StartNeighbor - StartOri
                    if Difference > 0 and Difference < NeighborDiff:
                        NeighborDiff = Difference
                        NeighborLocation = StartNeighbor
            elif OriginalGene.strand == "-":
                StartOri = OriginalGene.start
                for Gene in CandidateGenes:
                    StartNeighbor = Gene.start
                    Difference = StartOri - StartNeighbor
                    if Difference > 0 and Difference < NeighborDiff:
                        NeighborDiff = Difference
                        NeighborLocation = StartNeighbor
             
        else:
            raise Exception("direction must be defined as 'upstream' or 'downstream' with the quotes")
        NeighborGene=self.FindGeneByStart(NeighborLocation,Chromosome)
        return NeighborGene
    
    def FindIntergenicRegion(self, CNAG, strand):
        #self here is a list of all the genes as gene objects in Crypto genome
        OriginalGene = self.FindGeneByCNAG(CNAG)
        Chromosome = OriginalGene.chromosome
        CandidateGenes = self.FilterByChromosome(Chromosome)
        if strand == "+":
            StartIntergenic = 0
            EndIntergenic = OriginalGene.start
            for Gene in CandidateGenes:
                Start = Gene.end
                if int(Start) > int(StartIntergenic) and int(Start) < int(EndIntergenic):
                    StartIntergenic = Start
            #find upstream gene's end by going through all genes and finding which one is closest
        elif strand == "-":
            StartIntergenic=OriginalGene.end
            EndIntergenic = chromend[int(Chromosome) - 1]
            for Gene in CandidateGenes: 
                End = Gene.start
                if int(End) > int(StartIntergenic) and int(End) < int(EndIntergenic):
                    EndIntergenic = End           #find downstream gene's start by same method above 
        else:
            raise Exception("strand must be defined as '+' or '-' with the quotes")
        intergenic = Intergenic(CNAG, StartIntergenic, EndIntergenic, strand, Chromosome)
        return (intergenic)
         

class Gene(CommonEquality):
    """Represents genes by cnag and their corresponding position in crypto genome"""
    def __init__(self, CNAG, start, end, strand, chromosome):
        '''sets default values which don't make logical sense so I can tell if error happened'''
        self.CNAG = CNAG
        self.start = start
        self.end = end
        self.strand = strand
        self.chromosome = chromosome
    
    def __str__(self):
        '''for printing object as string'''
        return '%s %s %s %s %s' % (self.CNAG, self.start, self.end, self.strand, self.chromosome)
    
    def __len__(self):
        return abs(self.end-self.start)+1
 
class Peak(CommonEquality):
    '''Represents peaks by their name, position in crypto genome, and likelihood scores'''
    def __init__(self, name, chromosome, length, abs_summit_position, pileup, plogvalue, fold_enrichment, qlogvalue, distance_to_start_codon=None, matched=False):
        self.name = name
        self.chromosome = chromosome
        self.length = length
        self.abs_summit_position = abs_summit_position
        self.pileup = pileup
        self.plogvalue = plogvalue
        self.fold_enrichment = fold_enrichment
        self.qlogvalue = qlogvalue
        self.distance_to_start_codon= distance_to_start_codon
        self.matched = matched
        
    def __str__(self):
        return '%s %s %s %s %s %s %s %s' % (self.name, self.chromosome, self.abs_summit_position, self.pileup, self.plogvalue, self.fold_enrichment, self.qlogvalue, self.distance_to_start_codon)
        
    def FindDistanceToStartCodon(self, intergenic):
        if intergenic.strand == "+":
            self.distance_to_start_codon = int(intergenic.end) - int(self.abs_summit_position)
        elif intergenic.strand == "-":
            self.distance_to_start_codon = int(self.abs_summit_position) - int(intergenic.start)
        else: 
            print "Your intergenic object did not have a + or - strand, it is %s instead" % intergenic.strand

class Intergenic(CommonEquality):
    """Represents intergenic regions by the cnags for which they are promoters and their corresponding position in crypto genome"""
    def __init__(self, CNAG, start, end, strand, chromosome, type=None):
        self.CNAG = CNAG
        self.start = start
        self.end = end
        self.strand =  strand
        self.chromosome = chromosome
        self.length = int(end) - int(start)
        self.type = type
    
    def __str__(self):
        '''for printing object as string'''
        return '%s %s %s %s %s %s' % (self.CNAG, self.start, self.end, self.strand, self.chromosome, self.type)   

    def IsInIntergenic(self, peak):
        return int(self.start) < int(peak.abs_summit_position) and int(self.end) > int(peak.abs_summit_position) and self.chromosome == peak.chromosome
        

class PeaksList2(CommonEquality):
    
    def __init__(self, Peaks2=[]):
        self.Peaks2 = Peaks2
        
    def __len__(self):
        return len(self.Peaks2)
    
    def append(self, peak):
        self.Peaks2.append(peak)
        
    def __iter__(self):
        return iter(self.Peaks2)

class PeaksList(CommonEquality):
    
    def __init__(self, Peaks=[]):
     self.Peaks = Peaks

    def __iter__(self):
        return iter(self.Peaks)
    
    def __len__(self):
        return len(self.Peaks)
    
    def append(self, peak):
        self.Peaks.append(peak)
        
    def PopulateFromFile(self, filename):    
        with open(filename, 'r') as f:
            f.next()                                            #to skip first line of file which is a header
            for line in f:
                data=line.split('\t')
                chromname=data[Peak_chromosome_index]
                chrom=chromname.split('.')
                peakchromosome=chrom[1]    
                peaktemp = Peak(data[peak_name_index], peakchromosome, data[Peak_length_index], data[Abs_summit_index], data[Pileup_index], data[Log_p_value_index], data[Fold_enrichment_index], data[Log_q_value_index])  
                self.Peaks.append(peaktemp)

class IntergenicsList(CommonEquality):
    def __init__(self, Intergenics=[]):
        self.Intergenics = Intergenics
        
    def __iter__(self):
        return iter(self.Intergenics)
