# Scripts
## RNAseq_tools.py
Starting with fastq.gz files, aligns to genome with Tophat (using Bowtie 1), generates sorted, indexed bam files and counts reads in annotated transcripts using htseq-count. Generates all files necessary to proceed with DESeq2 analysis.
### Usage:  
```python RNAseq_tools.py directory number_of_threads <organism> <SE/PE>```
  
Required arguments:  
>directory : directory containing fastq.gz files (e.g. ./ or /data/jordan/JEB001/)  
>number_of_threads : processors to be used when running Tophat  
  
Optional arguments:  
>organism : default 'crypto' (C. neoformans H99), other options - 'pombe' (S. pombe) and 'cerevisiae' (S. cerevisiae S288C)  
>SE/PE : single end or paired end data, default - single end (SE)  
  
## ChIP_tools.py
Starting with fastq.gz files, trims adaptor with Cutadapt, aligns to genome with Bowtie 1 and generates sorted, indexed bam files. See log files for Cutadapt and Bowtie summaries.
### Usage:
```python ChIP_tools.py --directory fastq_directory --threads num_threads --organism crypto/pombe/cerevisiae/candida <--adaptor GATCGGAAGA> <--gff3 gff3_file> <--index bowtie_index_prefix>```  
  
Required positional arguments:  
>--directory : directory with fastq.gz files (e.g. ./ or /data/jordan/JEB001/)  
>--threads : number of processors to use  
>--organism : yeast genus or species - currently can be crypto/pombe/cerevisiae/candida  

Optional arguments:  
>--adaptor : change only if library preparation differs from typical Madhani lab ChIP-seq protocol  
>--gff3 and --index : can provide custom gff3 file and bowtie index i.e. for a different organism or updated versions of the genome. Must provide both.  
