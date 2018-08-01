# Introduction to GNU Make and Build Managers 
# Joanna Malukiewicz
# AGAR 2018 

## 1. Uses for GNU Make

Although Make is typically used to compile programs written in a language like C, it is also very handy for bioinformatics, data management and data analysis!

In bioinformatics and genomics, a typical pipeline could include steps to pre-treat and filter your sequencer files, perform _de novo_ or referenced genome alignment, further filtering, merging, variant calling, filtering, and then some type of analyses. 

You could write a series of shell scripts (which is fine), but Make and other build managers such as Snakemake and Scons give you much finer grain control and organization over running your pipeline than straight shell scripting. These build managers will let you 

  * Create a Makefile that describes how files in your pipeline depend on each other, and how to update out-of-date files.

  * Give flow between file dependencies and steps within your pipeline. 

  * Reduce redundancy to streamline your coding  (which will make your life much easier in the long run).

## 2. Example 1 (Simple)

  * Task:
    * Align a number of fastq Illumina files to the human reference genome with BWA 
    
```
$bwa mem sample1.F1.fastq sample2.R1.fastq > sample1.sam
```  
What is you had 10s or 100s files to process? This is where GNu Make can help. 
Start by setting up a _rule_ to process your task:

```
sample1.sam: sample1.F1.fastq sample2.R1.fastq
	bwa mem sample1.F1.fastq sample2.R1.fastq > sample1.sam
```

sample1.sam is the _target_ (or end result) of your rule 
the creation of the target is dependent on your input files, sample1.F1.fastq and sample2.R1.fastq, or your _dependencies_ . 
The second line, which always has to be indented, is called the _action line_ and contains the actual line of command. 

How do you modify the target rule to accomodate several fastq files? 

First we need to define a few variables 

```
# make a list of all the sam file targets you want to produce and define it as a variable 
DATABASE = output_files.txt 

# define the path to your reference genome as a variable 
REF= /path/to/ref/genome/hg38.fa

#define a variable that will act as a shell command that will read your list of output files via the shell cat command 
OUTPUTFILES := $(shell cat $(DATABASE))

```
	
Then set targets so that Make will process your entire list of output targets. 

To get Make to build everything at once, we have to introduce a phony target, which doesn't correspond to any files. 

Since a phony target doesn't actually exist, it won't ever be updated

```
all: $(OUTPUTFILES)

```
$(VARIALBLE_NAME) is a way to call whatever variable you set within a Makefile 

Then finally the generalized rule will use some automatic Make variables. 

```
    %.bam  :  %.unassembled.forward.fastq /%.unassembled.reverse.fastq
		bwa mem -t 10 $(REF) $^ > $@ 

```
% is a wild card character that will stand in for sample1, sample2, etc for the bam and fastq files 
$^ refers to your dependencies (%.unassembled.forward.fastq /%.unassembled.reverse.fastq)
$@ refers to your target (%.bam)

To process the Makefile you will type

```
$make all
```

All together, Make will go through your list of output bams and process them all one by one (though it it possible to parallelize the process).
If you type make all again, Make won't run again if it has processed all targets and rules. However, if you add in a new set of fastq files and run make, Make will just process the newly added fastq files. 

If you have multiple target rules, and you only want to run a rule for a specific target, write

```
$make <target_name>
```


The complete example code 
```
# make a list of all the sam file targets you want to produce and define it as a variable 
DATABASE = output_files.txt 

# define the path to your reference genome as a variable 
REF= /path/to/ref/genome/hg38.fa

#define a variable that will act as a shell command that will read your list of output files via the shell cat command 
OUTPUTFILES := $(shell cat $(DATABASE))

all: $(OUTPUTFILES)

%.sam  :  %.unassembled.forward.fastq /%.unassembled.reverse.fastq
	bwa mem -t 10 $(REF) $^ > $@  
```




## 3. Example 2  (More Advanced)

  * Task:
    * Download a series of cram files for the EPAS1 gene (chr2:46293667-46386703) from the 1000 Genomes database from a number of different individuals 
    * These cram files are already aligned to the human reference genome. 
    * We will merge these files and they can be used for functional or population genomics analyses. 
    
The command line task to obtain a cram file for a single sample is: 

```
samtools view -C -h ftp:/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CHS/HG00419/high_cov_alignment/HG00419.alt_bwamem_GRCh38DH.20150917.CHS.high_coverage.cram chr2:46293667-46386703 > EPAS1_indiv1.cram
```
Samtools is software for manipulating alignments of short genomic reads to a reference genome as either sam/bam/cram. Cram is a very compressed version of the file and -C indicates our preference for this format. -h indicates to include the "header" which is a lot of metadata before the actual alignments that may include the sequencing center, parts of the reference genome, samples, etc 

The command line for several samples is 

```
samtools view -C -h ftp:/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CHS/HG00419/high_cov_alignment/HG00419.alt_bwamem_GRCh38DH.20150917.CHS.high_coverage.cram chr2:46293667-46386703 > EPAS1_indiv1.cram
samtools view -C -h ftp:/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/GIH/NA20845/high_cov_alignment/NA20845.alt_bwamem_GRCh38DH.20150917.GIH.high_coverage.cram chr2:46293667-46386703> EPAS1_indiv2.cram
samtools view -C -h ftp:/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/ASW/NA19625/high_cov_alignment/NA19625.alt_bwamem_GRCh38DH.20150917.ASW.high_coverage.cram chr2:46293667-46386703 > EPAS1_indiv3.cram

```
What is you had 10s or 100s files to process? This is where GNu Make can help. 

```
#list of 1000 Genome individuals from which we will grab data 
input = /Users/arcova/Documents/Service/Bioinformatics_AAAGs/Makefiles/sample_list_short.txt #the := will expand only once whereas = will 

#define region of interest for EPAS1 gene, location obtained from ensembl.org for GRCh38.p12 human reference genome 
REGION= chr2:46293667-46386703 

#list of output crams 
output := /Users/arcova/Documents/Service/Bioinformatics_AAAGs/Makefiles/output_files_short.txt
OUTPUTFILES:= $(shell cat $(output))

.PHONY: all

all: $(OUTPUTFILES) rename

#get EPAS1 sequences from 1000 Genome dataset 
%.cram : 
	@echo '****' Extracting reads for $@ ...
	samtools view -C -h `grep $@ $(input)` $(REGION) > $@
	@echo '****' Success

#rename cram files to reflect that it is part of the EPAS1 gene 
rename:
	for file in *.cram ; do \
	mv $file ${file//high_coverage.cram/EPAS1.cram} \
	done

```


## Resources and References for this Presentation 


https://v4.software-carpentry.org/make/index.html

Greg Wilson: "Software Carpentry: Getting Scientists to Write Better
Code by Making Them More Productive".  Computing in Science &
Engineering, Nov-Dec 2006.

https://www.gnu.org/software/make/



