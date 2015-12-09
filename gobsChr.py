
#################################
#                               #
#   gobs module dependencies    #
#                               #
#################################

#import primer3
import HTSeq
import numpy
import copy
import pandas as pd
import types
import collections
import matplotlib.pyplot as plt
from copy import deepcopy



#################################
#                               #
#      gobs class definitions   #
#                               #
#################################

class GenomeSequence(object):
    def __init__(self,supercontigs):
        self.supercontigs=dict(enumerate(supercontigs))
        self.sizeDict={}
        for k,v in self.supercontigs.iteritems():
            self.sizeDict[str(k+1)]=len(v)
        print self.sizeDict
    def __len__(self):
        return len(self.supercontigs)
    def __getitem__(self,_i):
        return self.supercontigs[_i]
    def getSupercontig(self,chromosome=None):
        return self.supercontigs[chromosome]
    def getSequence(self,coord):
        return getsequence(self.supercontigs[coord.getChromosome()-1],coord.getStart(),coord.getEnd(),coord.getStrand())

class Coordinates(object):
    def __init__(self,chromosome,start,end,strand="Watson"):
        self.chromosome=chromosome
        self.start=start
        self.end=end
        self.strand=strand
        if self.start>self.end:
            print self.start,self.end
            raise ValueError("Illegal Coordinates")
    def __len__(self):
        return self.end-self.start+1
    def __repr__(self):
        return " ".join(str(x) for x in [self.chromosome,self.start,self.end, self.strand])
    def __lt__(self,other):
        if self.chromosome<other.chromosome: return True
        elif self.chromosome>other.chromosome: return False
        else:
            if self.start<other.start: return True
            elif self.start>=other.start: return False
            else: raise ValueError("Failed Comparison Error")
    def overlapsWith(self,other):
        if self.chromosome==other.chromosome:
            if self.start>other.start and self.start<other.end:
                return True
            elif self.end>other.start and self.end<other.end:
                return True
            elif other.start>self.start and other.start<self.end:
                return True
            elif other.end>self.start and other.end<self.end:
                return True
            else: return False
        else: return False
    def isWithin(self,other):
        if self.chromosome==other.chromosome:
            if self.start>other.start and self.end<other.end:
                return True
            else:
                return False
        else:
            return False
    def getChromosome(self):
        return self.chromosome
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end
    def getStrand(self):
        return self.strand
    def setStrand(self,new_strand):
        self.strand=new_strand

class Transcript(object):
    def __init__(self, transcript_name,chromosome,start_codon,stop_codon,strand,list_of_exons):
        self.name=transcript_name
        self.chromosome=chromosome
        self.start_codon=start_codon
        self.stop_codon=stop_codon
        self.strand=strand
        new_exons=[]
        for tpl in list_of_exons:
            if tpl[0]<tpl[1]:
                new_exons.append(tpl)
            else:
                new_exons.append((tpl[1],tpl[0]))
        sorted_exons=sorted(new_exons,key=lambda x: x[0])
        self.exons=sorted_exons
        self.five_prime_end=self.exons[0][0]
        self.three_prime_end=self.exons[-1][1]
        self.coord=Coordinates(self.chromosome,self.five_prime_end,self.three_prime_end,self.strand)
    def __repr__(self):
        return self.name
    def __len__(self):
        return abs(self.five_prime_end-self.three_prime_end)
    def getGeneName(self):
        return self.name[0:-2]
    def getName(self):
        return self.name
    def getCoordinates(self):
        return self.coord
    def getStrand(self):
        return self.strand
    def getCDSCoordinates(self):
        return Coordinates(self.chromosome,min(self.start_codon,self.stop_codon),max(self.start_codon,self.stop_codon),self.strand)
    def getNumberofIntrons(self):
        return len(self.exons)-1
    def getAllExonCoordinates(self):
        result=[]
        for exon in self.exons:
            start=min([exon[0],exon[1]])
            end=max([exon[0],exon[1]])
            result.append(Coordinates(self.chromosome, start,end,self.strand))
        return result
    def getAllIntronCoordinates(self):
        exons=self.getAllExonCoordinates()
        if self.strand=="Crick":
            exons=exons[::-1]
        introns=[]
        for (n, current_exon) in enumerate(exons):
            l=len(exons)-1
            if n<l:
                next_exon=exons[n+1]
                if self.strand=="Watson":
                    intron=Coordinates(self.chromosome,current_exon.getEnd()+1,next_exon.getStart()-1,self.strand)
                elif self.strand=="Crick":
                    intron=Coordinates(self.chromosome,next_exon.getEnd()+1,current_exon.getStart()-1,self.strand)
                introns.append(intron)
        return introns

    def getSplicedRNASequence(self,genome):
        result=''
        for exon in self.getAllExonCoordinates():
            result=result+genome.getSequence(exon)
            #print "Exon:", exon
        return result

    def getUnsplicedRNASequence(self,genome):
        self.Coordinates.getSequence(genome)

    def getIntronSequences(self,genome):
        sequence = []
        for intron in self.getAllIntronCoordinates():
            sequence.append(genome.getSequence(intron))
        return sequence

    def getSplicingAnalysisPrimerSet(self,genome,offset,primer_size,wigroom,opttemp):
        exons=self.getAllExonCoordinates()
        l=len(exons)
        primer_list=[]
        for n in range(len(exons)+1):
            if n<l-1:
                five_prime_primer=make_primer(genome,exons[n],offset,from_start=False,length=primer_size, wiggleroom=wigroom, optimalTm=opttemp)
                three_prime_primer=reversecomplement(make_primer(genome,exons[n+1],offset,from_start=True,length=primer_size,wiggleroom=wigroom,optimalTm=opttemp))
                if len(five_prime_primer)< primer_size-wigroom or len(three_prime_primer)<primer_size-wigroom:
                    primer_list.append(("None",0,"None",0))
                else:
                    exon1_sequence=genome.getSequence(exons[n])
                    exon2_sequence=genome.getSequence(exons[n+1])
                    five_prime_primer_location=exon1_sequence.find(five_prime_primer)
                    offset_of_five_prime_primer=len(exon1_sequence)-five_prime_primer_location
                    three_prime_primer_location=exon2_sequence.find(reversecomplement(three_prime_primer))+len(three_prime_primer)
                    offset_of_three_prime_primer=three_prime_primer_location
                    result=(five_prime_primer,offset_of_five_prime_primer,three_prime_primer, offset_of_three_prime_primer)
                    primer_list.append(result)
        return primer_list


    def getExonCoverage(self,cvg_genomic_array):
        translation={"1":"I","2":"II","3":"III"}
        result=0
        exons=[fixChromCoord(HTSeq_iv(x),translation) for x in self.getAllExonCoordinates()]
        print "Exons:",exons
        for exon in exons:
            result+=numpy.fromiter(cvg_genomic_array[exon],dtype='i').sum()
        return result

    def getIntronCoverage(self,cvg_genomic_array):
        translation={"1":"I","2":"II","3":"III"}
        result=0
        introns=[fixChromCoord(HTSeq_iv(x),translation) for x in self.getAllIntronCoordinates()]
        print "Introns:",introns
        if len(introns)<>0:
            for intron in introns:
                try:
                    result+=numpy.fromiter(cvg_genomic_array[intron],dtype='i').sum()
                except:
                    print "Could not get coverage for intron ",intron
                    return None
            return result
        else: return None


    def getFiveprimeUTRCoordinates(self):
        if self.strand=="Watson":
            print self.chromosome,self.getCoordinates().getStart(),self.start_codon,self.strand
            return Coordinates(self.chromosome,self.getCoordinates().getStart(),self.start_codon,self.strand)
        elif self.strand=="Crick":
            print self.name,self.chromosome,self.start_codon,self.getCoordinates().getEnd(),self.strand
            return Coordinates(self.chromosome,self.start_codon,self.getCoordinates().getEnd(),self.strand)
    def getThreeprimeUTRCoordinates(self):
        if self.strand=="Watson":
            return Coordinates(self.chromosome,self.stop_codon,self.getCoordinates().getEnd(),self.strand)
        elif self.strand=="Crick":
            return Coordinates(self.chromosome,self.getCoordinates().getStart(),self.stop_codon,self.strand)
    def getSplicedCDSRNACoordinates(self):
        pass
    def getProteinSequence(self):
        pass


class Gene(object):
    all_genes=set()
    def __init__(self,gene_name):
        self.name=gene_name
        Gene.all_genes.add(self.name)
        self.transcripts=[]
        self.upstream_element=None
        self.downstream_element=None
        self.coordinates=None
        self.strand=None
        self.literature_name=None
        self.annotation=None
        self.proteindomains=None
    def computeGeneBoundaryCoordinates(self):
        self.transcripts.sort(key=(lambda x: x.getCoordinates()))
        _all_exon_coordinates=[]
        for t in self.transcripts:
                _next_coords=t.getAllExonCoordinates()
                _all_exon_coordinates.extend(_next_coords)
        _first_transcript=self.transcripts[0]
        self.strand=_first_transcript.getCoordinates().getStrand()
        _all_exon_coordinates.sort()
        self.five_prime_most_end=_all_exon_coordinates[0].getStart()
        self.three_prime_most_end=_all_exon_coordinates[-1].getEnd()
        self.chromosome=_first_transcript.getCoordinates().getChromosome()
        self.coordinates=Coordinates(self.chromosome,self.five_prime_most_end,self.three_prime_most_end,self.strand)
        return self.coordinates
    def __repr__(self):
        return self.name
    def getCoordinates(self):
        return self.coordinates
    def getName(self):
        return self.name
    def addTranscript(self,transcript):
        self.transcripts.append(transcript)
    def getTranscripts(self):
        return self.transcripts
    def setUpstreamElement(self,upstream_element):
        self.upstream_element=upstream_element
    def setDownstreamElement(self,downstream_element):
        self.downstream_element=downstream_element
    def getUpstreamElement(self):
        return self.upstream_element
    def getDownstreamElement(self):
        return self.downstream_element
    def getPromoter(self):
        if self.strand=="Watson":
            coords=self.upstream_element.getCoordinates()
            self.promoter=IntergenicElement(coords)
        else:
            coords=copy.deepcopy(self.downstream_element.getCoordinates())
            coords.setStrand("Crick")
            self.promoter=IntergenicElement(coords)
        self.promoter.setName("p"+self.name)
        return self.promoter
    def getTerminator(self):
        if self.strand=="Watson":
            coords=self.downstream_element.getCoordinates()
            self.terminator=IntergenicElement(coords)
        else:
            coords=copy.deepcopy(self.upstream_element.getCoordinates())
            coords.setStrand("Crick")
            self.terminator=IntergenicElement(coords)
        self.terminator.setName("t"+self.name)
        return self.terminator
    def setAnnotation(self,annotation):
        self.annotation=annotation
    def setLiteratureName (self,literaturename):
        self.literature_name=literaturename
    def setProteinDomains(self,proteindomains):
        self.proteindomains=proteindomains
    def getAnnotation(self):
        return self.annotation
    def getLiteratureName(self):
        return self.literature_name
    def getProteinDomains(self):
        return self.proteindomains
    def getChromosomeLength(self, genome_obj):
        return len(genome_obj[self.coordinates.getChromosome()-1])
    def getDistancetoTelomere(self,genome_obj):
        distance_to_right_end=self.getChromosomeLength(genome_obj)-self.getCoordinates().getStart()
        distance_to_left_end=self.getCoordinates().getStart()
        if distance_to_right_end<distance_to_left_end:
            return distance_to_right_end
        else:
            return distance_to_left_end

class IntergenicElement(object):
    def __init__(self,coordinates):
        self.coordinates=coordinates
        self.name="generic intergenic element"
    def setName(self,name):
        self.name=name
    def __repr__(self):
        return self.name
    def __len__(self):
        return len(self.coordinates)
    def getCoordinates(self):
        return self.coordinates

class Centromere(IntergenicElement):
    def __init__(self,coordinates):
        self.coordinates=coordinates
        self.name="CEN "+str(self.coordinates.getChromosome())
    def __repr__(self):
        return self.name

class Chromosome(object):
    def __init__(self,Supercontig_number,genome_sequence):
        self.Supercontig_number=Supercontig_number
        self.genome_sequence=genome_sequence
        self.elements=[]
    def __len__(self):
        return len(self.genome_sequence[self.Supercontig_number-1])
    def __repr__(self):
        return str(self.Supercontig_number)
    def getName(self):
        return self.Supercontig_number
    def addElement(self,element):
        self.elements.append(element)
    def sortElements(self):
        return sorted(self.elements, key=self.getCoordinates())
    def getElements(self):
        return self.elements



#################################
#                               #
#      gobs function library    #
#                               #
#################################

def fixChromCoord(HTSeq_coord_obj,di):
    new=copy.deepcopy(HTSeq_coord_obj)
    correction=di[HTSeq_coord_obj.chrom]
    new.chrom=correction
    return new

def strip_except(text,allowed):
    return "".join([x for x in text if x in allowed])

def readfasta(fasta_file, alphabet):
    #This function reads in a fasta file of any kind and returns a list of 2 - tuples containing the description line the sequence
    f=open(fasta_file,"r")
    contents=f.read()
    f.close()
    pieces=contents.split(">")
    pieces=pieces[1:]
    li=[]
    for item in pieces:
        x=item.split("\n")
        gene_name=x[0]
        other_stuff=x[1:]
        text_of_other_stuff="".join(other_stuff)
        cleaned_up_text=strip_except(text_of_other_stuff,alphabet)
        li.append((gene_name,cleaned_up_text))
    return li

def getsequence(seq,start,end,strand):
    if strand=="Watson":
        return seq[start-1:end]
    elif strand=="Crick":
        return reversecomplement(seq[start-1:end])

def reversecomplement(seq):
    seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
    return "".join([seq_dict[base] for base in reversed(seq)])

def getsupercontig(supercontiglist, chromosome_number):
    return supercontiglist[chromosome_number-1]

def buildsupercontiglist (srcfilepath):
    li=readfasta(srcfilepath, "ACTGN")
    return ([x[1] for x in li])

def buildpfamdic(srcfilename):
    di={}
    with open(srcfilename,'r') as f:
        for line in f:
            contents=line.split('\t')
            gene_name=contents[1]
            pfamid=contents[3]
            pfamabbreviation=contents[4]
            pfamannotation=contents[5]
            if pfamid.startswith("PF"):
                if gene_name in di.keys():
                    di[gene_name].append((pfamid,pfamabbreviation,pfamannotation))
                else:
                    di[gene_name]=[(pfamid,pfamabbreviation,pfamannotation)]
    return di


def buildPombetranscriptdic(srcfilepath):
    #This function takes as its input PombeBase 12/2014 GFF3 file as input
    with open(srcfilepath,'r') as f:
        blocks=f.read().split('\n\n') #assumes that blocks for transcripts are separated by double newlines
        di={}
        for block in blocks:
            lines=block.split("\n")
            if len(lines)<2: continue
            lines_containing_mRNA=filter(lambda x: x.count("transcript")>0 ,lines)
            lines_containing_CDS=filter(lambda x: x.count("CDS")>0,lines)
            if len(lines_containing_mRNA)==0 or len(lines_containing_CDS)==0: continue
            if lines_containing_mRNA[0].split("\t")[0] not in ["I","II","III"]:continue
            for line_containing_mRNA in lines_containing_mRNA:
                line_pieces=line_containing_mRNA.split("\t")
                ninth_column=line_pieces[8]
                ninth_column_pieces=ninth_column.split(";")
                ninth_column_di={x.split("=")[0]:x.split("=")[1] for x in ninth_column_pieces}
                transcript_name=ninth_column_di["PARENT"].split(":")[1]
                lines_containing_exons=filter(lambda y: y.count("exon")>0,lines)
                exons=[]
                for line_containing_exon in lines_containing_exons:
                    pieces=line_containing_exon.split("\t")
                    chromosome=pieces[0]
                    if chromosome=="I":  chromosome=1
                    elif chromosome=="II":chromosome=2
                    elif chromosome=="III":chromosome=3
                    else: continue
                    start=int(pieces[3])
                    end=int(pieces[4])
                    strand=pieces[6]
                    if strand=="+":strand="Watson"
                    if strand=="-":strand="Crick"
                    exons.append((chromosome,start,end,strand))
            CDS_li=[]
            for line_containing_CDS in lines_containing_CDS:
                pieces=line_containing_CDS.split("\t")
                start=int(pieces[3])
                end=int(pieces[4])
                CDS_li.append(start)
                CDS_li.append(end)
            CDS_li.sort()
            if strand=="Watson":
                start_codon_pos=CDS_li[0]
                stop_codon_pos=CDS_li[-1]
            elif strand=="Crick":
                start_codon_pos=CDS_li[-1]
                stop_codon_pos=CDS_li[0]
            CDS_coord=[chromosome,start_codon_pos,stop_codon_pos,strand]
            di[transcript_name]=(CDS_coord,exons)
    print "done"
    return di


def buildH99transcriptdic(srcfilepath):
    print "Parsing GFF3 file: ",srcfilepath
    #This function takes as its input Broad H99 GFF3 file as input with "chr1" type chromosome coordinates
    #The function first makes a transcript isoform dictionary where the keys are the transcript isoform names
    # and the values are a list of text lines from the GFF3 file containing that name
    # The function then iterates through the transcript isoform dictionary, collecting lines that contain
    # the 'exon' or 'CDS' terms in the third column of the GFF3 file and then using that information to produce
    # exon coordinates along with the location of the start and stop codon for that transcript
    f=open(srcfilepath,'r')
    isoform_di={}
    result_di={}
    lines_containing_exons=[]
    lines_containing_CDS=[]
    for line in f:
        pieces=line.split("\t")
        if len(pieces)<3:
            continue
        if pieces[2]<>"exon" and pieces[2]<>"CDS":continue
        ninth_column=pieces[8]
        ninth_column_pieces=ninth_column.split(";")
        ninth_column_di={x.split("=")[0]:x.split("=")[1] for x in ninth_column_pieces}
        try:
            name=ninth_column_di["Parent"].strip()
        except:
            continue #No "Parent" entry for this line so continue to next line in file
        if name in isoform_di.keys():
            isoform_di[name].append(line)
        else:
            isoform_di[name]=[line]
    print len(isoform_di)," transcripts found in GFF3 file"
    f.close()

    for name,lines in isoform_di.iteritems():
        exons=[]
        CDS_li=[]
        lines_containing_exons=filter(lambda x: x.count("exon")>0 ,lines)
        lines_containing_CDS=filter(lambda x: x.count("CDS")>0,lines)
        #Collect the exons
        for line_containing_exon in lines_containing_exons:
            pieces=line_containing_exon.split("\t")
            chromosome=int(pieces[0].split("Chr_")[1])
            start=int(pieces[3])
            end=int(pieces[4])
            strand=pieces[6]
            if strand=="+":strand="Watson"
            if strand=="-":strand="Crick"
            exons.append([chromosome,start,end,strand])
        #Identify the start and stop codon coordinates
        for line_containing_CDS in lines_containing_CDS:
            pieces=line_containing_CDS.split("\t")
            ninth_column=pieces[8]
            ninth_column_pieces=ninth_column.split(";")
            ninth_column_di={x.split("=")[0]:x.split("=")[1] for x in ninth_column_pieces}
            start=int(pieces[3])
            end=int(pieces[4])
            CDS_li.append(start)
            CDS_li.append(end)
        CDS_li.sort()
        if strand=="Watson":
            start_codon_pos=CDS_li[0]
            stop_codon_pos=CDS_li[-1]
        elif strand=="Crick":
            start_codon_pos=CDS_li[-1]
            stop_codon_pos=CDS_li[0]
        CDS_coord=[chromosome,start_codon_pos,stop_codon_pos,strand]
        result_di[name]=(CDS_coord,exons)
    print "done building the H99 transcript dictionary"
    return result_di


def make_primer(genome_seq_obj,coord,offset,from_start=True,length=25,wiggleroom=5,optimalTm=60):
    sequence=genome_seq_obj.getSequence(coord)
    primer_list=[]
    posrange=range(wiggleroom)
    negrange=[-x for x in posrange]
    therange=negrange+posrange
    for wiggle in therange:
        if from_start==True:
            slice=sequence[offset+wiggle:offset+length+wiggle]
        else:
            slice=sequence[-offset-length+wiggle:-offset+wiggle]
#        tm=primer3.calcTm(slice)
        primer_list.append((slice,abs(tm-optimalTm)))
    primer_list.sort(key=lambda x:x[1])
    best_primer=primer_list[0][0]
    primer_list=[]
    for n in range(wiggleroom):
        slice1=best_primer[:-n]
        slice2=best_primer[n:]
#        primer_list.append((slice1,abs(primer3.calcTm(slice1)-optimalTm)))
#        primer_list.append((slice2,abs(primer3.calcTm(slice2)-optimalTm)))
    primer_list.sort(key=lambda x:x[1])
    best_best_primer=primer_list[0][0]
    return best_best_primer

def generate_plate_labels():
    rows=["A","B","C","D","E","F","G","H"]
    cols=[str(x+1) for x in range(12)]
    while True:
        for r in rows:
            for c in cols:
                yield r+c

def get_configuration(config_file):
    di={}
    with open(config_file,'r') as f:
        for line in f:
            k,v= [x.strip() for x in line.split("=")]
            di[k]=v
    return di


def purge(setofgenomicpositions, setofgenomicfeatures,halfwinwidth):
    interval_dict={}
    for pos in setofgenomicpositions:
        interval_dict[pos]=HTSeq.GenomicInterval(pos.chrom,pos.start-halfwinwidth,pos.end+halfwinwidth,pos.strand)
    for feature in setofgenomicfeatures:
        for position,interval in interval_dict.iteritems():
            if interval.overlaps(feature.iv):
                try:
                    setofgenomicpositions.remove(position)
                except:
                    pass
    return True


def smooth(x,window_len=201,window='flat'):
    """smooth the data using a window with requested size.
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    example:
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    See also:
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=numpy.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')
    y=numpy.convolve(w/w.sum(),s,mode='valid')
    return y

def ceilingize(vector,ceiling):
    result=numpy.zeros(len(vector))
    for value in range(len(vector)):
        if vector[value]>ceiling:
            result[value]=ceiling
        else:
            result[value]=vector[value]
    return result

def buildCoverageVector((bam_file_object,genome_size_di)):
    coverage=HTSeq.GenomicArray(genome_size_di,stranded=True,storage="ndarray",typecode="i")
    print "Building coverage vector for ",bam_file_object.filename
    for almnt in bam_file_object:
        if almnt.aligned:
            if almnt.iv.chrom=="I" or almnt.iv.chrom=="II" or almnt.iv.chrom=="III":
                if coverage[almnt.iv]:
                    coverage[almnt.iv] +=1
                else:
                    print "no coverage data for ",almnt.iv
    return coverage


def buildFeatureProfile(set_of_positions,cvg_vector,ceiling,halfwinwidth):
    sense_profile =numpy.zeros(2*halfwinwidth,dtype="i")
    antisense_profile=numpy.zeros(2*halfwinwidth,dtype="i")
    for p in set_of_positions:
        sense_window=HTSeq.GenomicInterval(p.chrom,p.pos-halfwinwidth,p.pos+halfwinwidth,"+")
        antisense_window=HTSeq.GenomicInterval(p.chrom,p.pos-halfwinwidth,p.pos+halfwinwidth,"-")
        try:
            sense_wincvg=numpy.fromiter(cvg_vector[sense_window],dtype="i",count=2*halfwinwidth)
            antisense_wincvg=numpy.fromiter(cvg_vector[antisense_window],dtype="i",count=2*halfwinwidth)
        except:
            print "had a problem with feature ",p
            continue
        if ceiling<>False:
            sense_wincvg=ceilingize(sense_wincvg,ceiling)
            antisense_wincvg=ceilingize(antisense_wincvg,ceiling)
        sense_profile += sense_wincvg
        antisense_profile += antisense_wincvg
        antisense_profile=antisense_profile[::-1]
    return (sense_profile,antisense_profile)

def buildFeatureProfile2(set_of_positions,sortedbamfile,ceiling,halfwinwidth):
    sense_profile =numpy.zeros(2*halfwinwidth,dtype="i")
    antisense_profile=numpy.zeros(2*halfwinwidth,dtype="i")
    for p in set_of_positions:
        sense_cvg_result =numpy.zeros(2*halfwinwidth,dtype="i")
        antisense_cvg_result=numpy.zeros(2*halfwinwidth,dtype="i")
        window=HTSeq.GenomicInterval(p.chrom,p.pos-halfwinwidth,p.pos+halfwinwidth,".")
        try:
            for almnt in sortedbamfile[window]:
                if p.strand == "+":
                    start_in_window = almnt.iv.start - p.pos + halfwinwidth
                    end_in_window   = almnt.iv.end   - p.pos + halfwinwidth
                    sense_cvg_result[ start_in_window : end_in_window ] += 1
                else:
                    start_in_window = p.pos + halfwinwidth - almnt.iv.end
                    end_in_window   = p.pos + halfwinwidth - almnt.iv.start
                    antisense_cvg_result[ start_in_window : end_in_window ] += 1
        except:
            print "Had a problem with position ",p
        sense_cvg_result=ceilingize(sense_cvg_result,ceiling)
        antisense_cvg_result=ceilingize(antisense_cvg_result,ceiling)
        sense_profile=sense_profile+sense_cvg_result
        antisense_profile=antisense_profile+antisense_cvg_result
    return (sense_profile,antisense_profile)


def buildFeatureProfile3(genomic_feature_set,alignment_type,cvg_vector,flank_size,slice_size,alignment_points=None):
    #alignment can be center, left, right, or variable
    #if variable, the function excepts kward 'alignment_points' to be a list of the indices of the alignment points
    features=list(genomic_feature_set)
    features.sort(key=lambda x: x.iv.length)
    for gf in features:
        try:
            data=numpy.fromiter(cvg_vector[gf.iv],dtype=int)
        except:
            print "No coverage for",gf.iv
            continue
        #if numpy.sum(data)==0:
        #   features.remove(gf)
    for gf in features:
        c=gf.iv.chrom
        start=int(gf.iv.start)
        end=int(gf.iv.end)
        strand=gf.iv.strand
        gf.iv.start=start-flank_size
        gf.iv.end=end+flank_size
    max_interval=max(x.iv.length for x in features)  #have to also code sort alignment points
    num_rows=len(features)
    if alignment_type<>"variable":
        num_cols=max_interval+2*flank_size
    elif alignment_type=="variable":
        num_cols=2*max_interval
        mid_col=max_interval
    profile=numpy.full((num_rows,num_cols),-10,dtype=int)
    index_li=[]
    for row,gf in enumerate(features):
        try:
            data=numpy.fromiter(cvg_vector[gf.iv],dtype=int)
        except:
            print "could not get data for ",gf
            continue
        #if gf.iv.strand=="-":
        #    data=data[::-1]
        print gf.name,gf.attr,gf.iv.length
        if type(gf.attr)==types.DictType:
            attribute=gf.attr["description"]
        else:
            attribute=gf.attr
        index_li.append(gf.name+":"+attribute+":"+str(gf.iv.length))
            #print "Could not get data for interval: ",gf
            #continue
        if alignment_type=="left": col_slice=slice(0,len(data))
        elif alignment_type=="right": col_slice=slice(num_cols-len(data),num_cols)
        elif alignment_type=="center": col_slice=slice(int(num_cols-len(data))/2,int(len(data)+((num_cols-len(data))/2)))
        elif alignment_type=="variable":
            alignment_point=alignment_points[row]
            col_slice=(mid_col-alignment_point,len(data)+(mid_col-alignment_point))
        profile[row,col_slice]=ceilingize(data,1)
    profile=profile[0:len(index_li),:]
    if alignment_type=="left":  profile=profile[:,0:(slice_size+flank_size)]
    elif alignment_type=="right":  profile=profile[:,-(slice_size+flank_size):]
    df=pd.DataFrame(profile,index=index_li,dtype=int)
    df.index.name="UID"
    print profile.shape
    return df

def buildFeaturePositionSet(feature_type,feature_attr_dict,GFF_reader_object,which_end):
    result=set()
    for feature in GFF_reader_object:
        if feature.type==feature_type:
            flag=True
            for k,v in feature_attr_dict.iteritems():
                if feature.attr[k]<>v:
                    flag=False
                    break
            if flag==True:
                if which_end=="start":
                    result.add(feature.iv.start_d_as_pos)
                if which_end=="end":
                    result.add(feature.iv.end_d_as_pos)
    return result

def buildFeatureSet(feature_type,feature_attr_dict,GFF_reader_object):
    result=set()
    for feature in GFF_reader_object:
        if feature.type==feature_type:
            flag=True
            for k,v in feature_attr_dict.iteritems():
                if feature.attr[k]<>v:
                    flag=False
                    break
            if flag==True:
                result.add(feature)
    return result

def HTSeq_iv(coordinate_object):
    #converts gobs Coordinate object to HTSeq GenomicInterval object
    chrom="Chr_"+str(coordinate_object.getChromosome())
    start=coordinate_object.getStart()-1
    end=coordinate_object.getEnd()
    strand=coordinate_object.getStrand()
    if strand=="Watson":strand="+"
    elif strand=="Crick":strand="-"
    else:
        raise ValueError,"Strand was not Watson or Crick"
    return HTSeq.GenomicInterval(str(chrom),start,end,strand)

def Coordinate_from_HT_Seq_Interval(HTSeq_iv):
    return Coordinates(HTSeq_iv.chrom,HTSeq_iv.start+1,HTSeq_iv.end,HTSeq_iv.strand)

def interval_reverse_complement(input_interval):
    cp=deepcopy(input_interval)
    if input_interval.strand=="+": cp.strand="-"
    if input_interval.strand=="-": cp.strand="+"
    return cp

def buildGenomeDatabase(configuration_file):
    ####################################################################################################
    #Use custom parsers in genomefunctions to build dicts of CDSs, supercontigs, and transcript models #
    ####################################################################################################
    config=get_configuration(configuration_file)
    species=config["SPECIES"]
    supercontigs=buildsupercontiglist(config["FASTA"])
    print "Supercontig lengths:"
    print [len(x) for x in supercontigs]
    print "SPECIES= ",species
    if species=="SCHIZOSACCHAROMYCESPOMBE":
        transcriptdi=buildPombetranscriptdic(config["GFF3"])
    elif species=="CRYPTOCOCCUSNEOFORMANSVARGRUBII":
        transcriptdi=buildH99transcriptdic(config["GFF3"])
    else:
        raise ValueError ("no species value given")
    ##########################
    #Build Centromere Objects#
    ##########################
    centromeres=[]
    with open(config["CEN"]) as f:
        print "Centromere coordinates: "
        for line in f:
            pieces=line.split(":")
            chromosome=int(pieces[0])
            start=int(pieces[1].split("-")[0])
            end=int(pieces[1].split("-")[1].strip())
            print chromosome,start,end
            centromeres.append(Centromere(Coordinates(chromosome,start,end,"Watson")))
    ###################################
    # Populate a GenomeSequence object#
    ###################################
    print "Building GenomeSequence object"
    genome=GenomeSequence(supercontigs)
    print len(genome)," Supercontigs loaded"
    ####################################
    #Build a List of Transcript Objects#
    ####################################
    print "Building Transcript objects"
    transcript_collection=[]
    for transcript,info in transcriptdi.iteritems():
        CDS_info=info[0]
        chromosome=CDS_info[0]
        start=CDS_info[1]
        stop=CDS_info[2]
        strand=CDS_info[3]
        exon_list=info[1]
        exon_starts_stops = []
        for x in exon_list:
            exon_starts_stops.append((int(x[1]), int(x[2])))
        transcript_obj=Transcript(transcript,chromosome,start,stop,strand,exon_starts_stops)
        transcript_collection.append(transcript_obj)
    ##################################
    #Test out some Transcript methods#
    ##################################
    print "Testing transcript methods"
    test_list=transcript_collection[0:10]
    for test_transcript in test_list:
        print test_transcript,"=TRANSCRIPT NAME"
        print test_transcript.getCoordinates(),"=TRANSCRIPT COORDINATES"
        print test_transcript.getCoordinates().getChromosome(),"= CHROMOSOME"
        print test_transcript.getCoordinates().getStart(),"=START"
        print test_transcript.getCoordinates().getEnd(),"=END"
        print test_transcript.getCoordinates().getStrand(),"=STRAND"
        print "ALL EXON COORDINATES :"
        print test_transcript.getAllExonCoordinates()
        print "EXON SEQUENCES:"
        print [genome.getSequence(x) for x in test_transcript.getAllExonCoordinates()]
        print "INTRON SEQUENCES:"
        intron_list=[genome.getSequence(x) for x in test_transcript.getAllIntronCoordinates()]
        print intron_list
        print [len(y) for y in intron_list ]
        print "INTRON LENGTHS"
        print [len(z) for z in test_transcript.getAllIntronCoordinates()]
        print "CDS SEQUENCES"
        print genome.getSequence(test_transcript.getCDSCoordinates())
        print "--------------------"


    ###############################
    #Build a list of Gene objects #
    ###############################
    print "Building Gene objects"
    gene_coll=[]
    for transcript in transcript_collection:
        transcript_gene_name=transcript.getGeneName()
        if not (transcript_gene_name in Gene.all_genes):
            gene_coll.append(Gene(transcript_gene_name))
        for gene in gene_coll:
            if gene.getName()==transcript_gene_name:
                gene.addTranscript(transcript)
                break
    for gene in gene_coll:
        gene.computeGeneBoundaryCoordinates()
    ############################
    #Test out some Gene methods#
    ############################
    print "Testing gene methods"
    gene_coll.sort(key=lambda x:x.getCoordinates())
    print [x.computeGeneBoundaryCoordinates() for x in gene_coll][0:5]
    #################################
    #Populate Chromosomes with Genes#
    #################################
    chromosomes=[]
    for num in range(len(genome)):
        chromosomes.append(Chromosome(num+1,genome))
    for gene in gene_coll:
        for chromosome in chromosomes:
            if gene.getCoordinates().getChromosome()==chromosome.getName():
                chromosome.addElement(gene)
    ########################################################################################
    #Build IntergenicElement Objects and set Upstream and Downstream Elements for each Gene#
    ########################################################################################
    print "Adding IntergenicElement objects to each GeneObject"
    intergenic_floor=100
    #if an intergenic region is less then a predefined value (including overlapping genes), set the floor to this value
    overlapcounter=0
    three_prime__to_three_prime_overlapcounter=0
    for chromosome in chromosomes:
        print "Adding intergenic elements to Chromosome: ",chromosome
        for n,gene in enumerate(chromosome.getElements()):
            if n==0:
                #first gene on a chromosome
                print "The first gene in chromosome ",chromosome," is ",gene.getName()
                print "and it's coordinates are: ", gene.getCoordinates()
                try:
                    gene.setUpstreamElement(IntergenicElement(Coordinates(chromosome.getName(),1,gene.getCoordinates().getStart()-1,"Watson")))
                except:
                    gene.setUpstreamElement(IntergenicElement(Coordinates(chromosome.getName(),1,2,"Watson")))
                    print "First gene runs to end of chromosome -- likely truncated"
            current_gene_end=gene.getCoordinates().getEnd()
            try:
                next_gene=chromosome.getElements()[n+1]
                next_gene_coordinates=next_gene.getCoordinates()
            except:
                #last gene on a chromosome
                print "last gene on chromosome ",chromosome," is ",gene.getName()
                try:
                    gene.setDownstreamElement(IntergenicElement(Coordinates(chromosome.getName(),gene.getCoordinates().getEnd()+1,len(chromosome),"Watson" )))
                except:
                    gene.setDownstreamElement(IntergenicElement(Coordinates(chromosome.getName(),len(chromosome)-1,len(chromosome),"Watson" )))
                    print "Last gene does not have a downstream intergenic region.  Likely runs off the end of the chromosome"
                continue #break out of current for loop iteration
            next_gene_start=next_gene_coordinates.getStart()
            #Deal with the special case of genes that overalp or are below a floor we set for a minimal promoter/terminator
            if next_gene_start<=current_gene_end+intergenic_floor:
            #print gene," and ",next_gene, "overlap.  Here are their coordinates: ",gene.getCoordinates()," vs ",next_gene_coordinates
                if gene.getCoordinates().getStrand()=="Watson" and next_gene_coordinates.getStrand()=="Crick":
                    three_prime__to_three_prime_overlapcounter=three_prime__to_three_prime_overlapcounter+1
                overlapcounter=overlapcounter+1
                downstream_region_coords=Coordinates(gene.getCoordinates().getChromosome(),current_gene_end+1, current_gene_end+intergenic_floor+1,"Watson")
                gene.setDownstreamElement(IntergenicElement(downstream_region_coords))
                next_gene.setUpstreamElement(IntergenicElement(downstream_region_coords))
                continue
            next_intergenic_coords=Coordinates(gene.getCoordinates().getChromosome(),current_gene_end+1,next_gene_start-1,"Watson")
            next_intergenic_element=IntergenicElement(next_intergenic_coords)
            gene.setDownstreamElement(next_intergenic_element)
            next_gene.setUpstreamElement(next_intergenic_element)
    print overlapcounter," cases of transcript overlap"
    print three_prime__to_three_prime_overlapcounter," cases of 3' to 3' overlap"
    #############################################################
    #Test out Gene methods that obtain promoters and terminators#
    #############################################################
    print "Testing Gene and Chromosome methods"
    print "There are ",len(chromosomes)," chromosomes"
    a,b,c,d=0,0,0,0
    for chromosome in chromosomes:
        list_of_genes=chromosome.getElements()
        list_of_promoters=[x.getPromoter() for x in list_of_genes]
        list_of_terminators=[y.getTerminator() for y in list_of_genes]
        longest_promoter=max(list_of_promoters, key=lambda x:len(x))
        longest_terminator=max(list_of_terminators, key=lambda x:len(x))
        print "For chromosome ",chromosome
        print "number of genes = ", len(chromosome.getElements())
        print "Longest Promoter = ", longest_promoter," ",len(longest_promoter)
        print "Longest Terminator =",longest_terminator," ", len(longest_terminator)
        for g in list_of_genes:
            for cen in centromeres:
                if g.getCoordinates().overlapsWith(cen.getCoordinates()):
                    print g,' overlaps with ',cen
                if g.getCoordinates().isWithin(cen.getCoordinates()):
                    print g,' is within ',cen
            if g.getDistancetoTelomere(genome)<40000:
                pass
                #print g,' is ',g.getDistancetoTelomere(genome),'bp from a telomere'
        for n, g in enumerate(list_of_genes[:-1]):
            strand1=g.getCoordinates().getStrand()
            next_gene=list_of_genes[n+1]
            strand2=next_gene.getCoordinates().getStrand()
            if strand1=="Watson" and strand2=="Watson": a+=1
            if strand1=="Crick" and strand2=="Crick":b+=1
            if strand1=="Crick" and strand2=="Watson":c+=1
            if strand1=="Watson" and strand2=="Crick": d+=1
    print ">>:",a
    print "<<:",b
    print "<>:",c
    print "><:",d
    print "total=",a+b+c+d

    return genome,chromosomes


def count_reads((bamfilereader,feature_di)):
    ##############################################################################################################################
    # This function returns a counter object with read counts for each item in the feature dictionary                            #
    # Inputs:  sortedbamfile is sorted .bam file                                                                                 #
    #          feature_di is a dictionary in which the keys are string names and the values are a HT-seq Genomic interval object  #
    # Dependencies:  HTSeq, collections modules                                                                                  #
    # Output:  a collections.Counter object which the keys are the keys of the feature dict and the values are counts            #
    ##############################################################################################################################
    print "Counting Reads in ",bamfilereader.filename,"..."
    cntr=0
    read_counter=collections.Counter()
    for name,interval in feature_di.iteritems():
        cntr=cntr+1
        if cntr%1000==0: print cntr, "features analyzed"
        for almt in bamfilereader[interval]:
            read_counter[name]+=1
    print sum(v for v in read_counter.values()), "Reads Counted"
    return read_counter

def calculate_read_position(feature,read, exons_only=True,introns=[]):
    #'feature' and 'read' are gobs Coordinate objects and read should be 'in' feature
    # returns distance from the start of feature to start of read
    # if exons_only flag is true then the length of introns in that interval are subtracted
    if read.getStrand()==feature.getStrand() and read.getChromosome()==feature.getChromosome() and read.isWithin(feature):
        if read.getStrand()=="Watson":
            ruler=Coordinates(read.getChromosome(), feature.getStart(),read.getStart(),strand="Watson")
        else:
            ruler=Coordinates(read.getChromosome(), read.getEnd(),feature.getEnd()+2,strand="Crick")
        introns_in_ruler=[]
        if exons_only==True and introns<>[]:
            for intron in introns:
                if intron.isWithin(ruler):
                    introns_in_ruler.append(intron)
            length_of_introns_in_ruler=0
            for item in introns_in_ruler:
                length_of_introns_in_ruler+=(len(item))
            return len(ruler)-length_of_introns_in_ruler
        else:
            return len(ruler)
    else:
        #print "strand and/or chromosome of read and feature do not match or read is not within feature"
        return -1

def get_alignments_in_an_interval(sortedbamfile,interval):
    #'interval' is a HTSeq 'GenomicInterval object
    #the function returns a 'list' of HTSeq 'Alignment' objects
    bamfile=HTSeq.BAM_Reader(sortedbamfile)
    return bamfile[interval]

def transcript_factory(config_file_path):
    genome,chromosomes=buildGenomeDatabase(config_file_path)
    for chromosome in chromosomes:
        for gene in chromosome.getElements():
            for transcript in gene.getTranscripts():
                yield transcript, genome

def simpleIOoptions():
    from optparse import OptionParser
    parser=OptionParser()
    parser.add_option("-i")
    parser.add_option("-o")
    (opts,args)=parser.parse_args()
    print opts
    inputfile=opts.i
    outputfileprefix=opts.o
    return inputfile,outputfileprefix

def barplot(data, x_axis_label,y_axis_label):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    rects=ax.bar(range(len(data)),data,color='blue')
    ax.set_xlabel(x_axis_label)
    ax.set_ylabel(y_axis_label)
    return ax

def IOoptions(list_of_options):
    from optparse import OptionParser
    parser=OptionParser()
    for item in list_of_options:
        parser.add_option(item)
    (opts,args)=parser.parse_args()
    return opts

def count_reads_in_a_interval(bamfilereader,interval):
    cntr=0
    for almt in bamfilereader[interval]:
        cntr=cntr+1
    return cntr

def count_reads__whose_start_matches_the_start_of_an_interval(bamfilereader,exon_interval_rc):
    cntr=0
    interval=HTSeq.GenomicInterval(chrom=exon_interval_rc.chrom,start=exon_interval_rc.start-50,end=exon_interval_rc.end+50,strand=exon_interval_rc.strand)
    try:
        for almt in bamfilereader[interval]:
            if almt.iv.strand==exon_interval_rc.strand:
                if exon_interval_rc.strand=="+":
                    if almt.iv.start-exon_interval_rc.start==0:
                        cntr=cntr+1
                elif exon_interval_rc.strand=="-":
                    if almt.iv.end-exon_interval_rc.end==0:
                        cntr=cntr+1
    except: print "Failed count for exon interval reverse complement: ", exon_interval_rc
    return cntr

def count_splicing_intermediates((bamfilereader,exon_di)):
    read_counter=collections.Counter()
    for name,exons in exon_di.iteritems():
        for num,exon_iv in enumerate(exons[0:-1]):
            exon_interval_rc=interval_reverse_complement(exon_iv)
            read_counter[(name,num+1)]=count_reads__whose_start_matches_the_start_of_an_interval(bamfilereader,exon_interval_rc)
    print len(read_counter)," exon three-prime ends analyzed"
    return read_counter
















