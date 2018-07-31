## Creating fastqs

Files prepared by Maria Nieves Colon. Commands based, in part, on tutorial by T. Webster available at https://github.com/thw17/ASU_BIO543_Genome_Assembly_2018/tree/master/fastq

These fastqs were created by subsetting data from 1000 genomes low coverage bams (Example: NA18501.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam) available from here: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/

We randomly selected six individuals from three 1000 genomes populations: CEU, YRI and PUR (N=18). To find the names of the individuals we made a popinfo file and searched it using grep and head:
```
grep "CEU" 1KGP3.popinfo.txt | head -n6
NA06984 NA06984 CEU EUR 1
NA06985 NA06985 CEU EUR 2
NA06986 NA06986 CEU EUR 1
NA06989 NA06989 CEU EUR 2
NA06994 NA06994 CEU EUR 1
NA07000 NA07000 CEU EUR 2

grep "YRI" 1KGP3.popinfo.txt | head -n6
NA18486 NA18486 YRI AFR 1
NA18488 NA18488 YRI AFR 2
NA18489 NA18489 YRI AFR 2
NA18498 NA18498 YRI AFR 1
NA18499 NA18499 YRI AFR 2
NA18501 NA18501 YRI AFR 1  

grep "PUR" 1KGP3.popinfo.txt | head -n6
HG00551 HG00551 PUR AMR 2
HG00553 HG00553 PUR AMR 1
HG00554 HG00554 PUR AMR 2
HG00637 HG00637 PUR AMR 1
HG00638 HG00638 PUR AMR 2
HG00640 HG00640 PUR AMR 1
```

The bam files were downloaded from the 1000 genomes ftp website using their download tool Aspera (see info on Aspera here: http://www.internationalgenome.org/category/aspera/). Example command to download .bam file, note that we also downloaded indexed .bai file for each sample:
```
./ascp -i asperaweb_id_dsa.openssh -Tr -Q -l 100M -P33001 -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:/vol1/ftp/phase3/data/NA18501/alignment/NA18501.mapped.ILLUMINA.bwa.YRI.low_coverage.20130415.bam .
```
To extract reads mapped to MT, we used the command:

```
samtools view -b NA18501.mapped.ILLUMINA.bwa.CEU.low_coverage.20120522.bam MT | samtools bam2fq -1 NA18501_R1.fq -2 NA18501_R2.fq -
```

These reads need to be sorted and read pairing needs to be restored, so we used repair.sh available from bbmap:
```
/Users/Maria/Install/bbmap/repair.sh in1=NA18501_R1.fq in2=NA18501_R2.fq out1=NA18501_R1.fastq out2=NA18501_R2.fastq
```

Finally, we subset a total of 10,000 paired reads from each sample with the following commands. We also renamed the fastq as POP_Sample.
The files were then compressed using gzip to save space:

```
seqtk sample -s100 NA18501_R1.fastq 5000 > YRI_NA18501_MT.R1.fastq
seqtk sample -s100 NA18501_R2.fastq 5000 > YRI_NA18501_MT.R2.fastq

gzip *.fastq
```
