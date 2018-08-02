#simple make file to retrive EPAS1 and turn it into its CDS 

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

