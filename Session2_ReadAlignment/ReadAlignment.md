
## Sessions 2 and 3: Assembly, variant calling, variant filtering

### File formats covered
fasta
fastq
sam/bam
vcf

### Step 1 -- What is a reference genome?
- what is a reference?
- what is fasta format?
- how do we prepare a reference for read mapping?
  - bwa index, fai index (samtools faidx), sequence dictionary


### Step 1 commands 

Note: ensure Miniconda or Anaconda is installed first. Also, the conda config commands *must be run in this order*

```
### Configure channel priorities for Conda
conda config --add channels r

conda config --add channels defaults

conda config --add channels conda-forge

conda config --add channels bioconda

### Create Conda environment for analyses
conda create -n agar2018 python=3.6 snakemake samtools bwa bioawk fastqc multiqc bbmap qualimap gatk4 vcftools picard

### Load Conda environment
source activate agar2018

### Index reference genome
samtools faidx reference/human_v37_MT.fasta

samtools dict -o reference/human_v37_MT.dict reference/human_v37_MT.fasta

bwa index reference/human_v37_MT.fasta

### Fastqc
# Example command for one file
fastqc -o fastqc_results fastq/CEU_NA07000_MT.R1.fastq.gz

### Multiqc
export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && multiqc --interactive -o multiqc_results fastqc_results"

### Fastq Trimming
# Example command for one sample
bbduk.sh -Xmx1g in1=fastq/PUR_HG00553_MT.R1.fastq.gz in2=fastq/PUR_HG00553_MT.R2.fastq.gz out1=trimmed_fastqs/PUR_HG00553_trimmed_read1.fastq.gz out2=trimmed_fastqs/PUR_HG00553_trimmed_read2.fastq.gz ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe qtrim=rl trimq=15 minlen=50 maq=20

### Trimmed Fastqc
# Example for one sample
fastqc -o fastqc_trimmed_results trimmed_fastqs/CEU_NA06986_trimmed_read1.fastq.gz trimmed_fastqs/CEU_NA06986_trimmed_read2.fastq.gz

### Trimmed Multiqc
export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && {params.multiqc} -o multiqc_trimmed_results fastqc_trimmed_results"
```

### Step 2 -- Fastq quality control
- what is a fastq?
- basics of exploring a fastq
  - bioawk
- how do we explore fastq quality?
  - fastqc
  - multiqc
- how do we trim fastqs?
  - bbduk (also trimmomatic/trim galore/etc.)
  - why and when do we trim fastqs?
    - what is an adapter?

```
### Mapping reads (also fixmate and sorting)
# Example for one sample
bwa mem -R '@RG\tID:YRI_NA18498\tSM:YRI_NA18498' reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz| samtools fixmate -O bam - - | samtools sort -O bam -o bams/YRI_NA18498.human_v37_MT.sorted.bam

# Index mapped & sorted BAM (one sample example)
samtools index bams/YRI_NA18498.human_v37_MT.sorted.bam

### Marking duplicates
# Example for one sample
picard -Xmx1g MarkDuplicates I=bams/PUR_HG00640.human_v37_MT.sorted.bam O=bams/PUR_HG00640.human_v37_MT.sorted.mkdup.bam M=stats/PUR_HG00640.human_v37_MT.picard_mkdup_metrics.txt

# Index mkdup BAM (one sample example)
samtools index bams/PUR_HG00640.human_v37_MT.sorted.mkdup.bam

### BAM statistics
# BAM stats (one sample example)
samtools stats bams/CEU_NA06994.human_v37_MT.sorted.mkdup.bam | grep ^SN | cut -f 2- > stats/CEU_NA06994.human_v37_MT.sorted.mkdup.bam.stats

# qualimap
# first, run per sample (one sample example)
qualimap bamqc -bam bams/CEU_NA06986.human_v37_MT.sorted.mkdup.bam -nt 1 -outdir stats/qualimap/human_v37_MT/CEU_NA06986/

# then run qualimap's multi-bamqc on all individual results
qualimap multi-bamqc -d stats/qualimap/human_v37_MT/qualimap.list -outdir stats/qualimap/human_v37_MT/
```








