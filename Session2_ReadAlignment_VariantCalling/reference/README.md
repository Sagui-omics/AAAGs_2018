### Reference information
Reference contains the mitochondrial (MT) sequence from the 1000 genomes reference (human_g1k_v37.fasta) available here: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/

The MT sequence was extracted with the command:
```
samtools faidx human_g1k_v37.fasta MT > human_v37_MT.fasta
```

### Commands to index reference
For our pipeline, you'll need a few index files that can be easily generated with the following three commands:
```
samtools faidx reference/human_v37_MT.fasta

samtools dict -o reference/human_v37_MT.dict reference/human_v37_MT.fasta

bwa index reference/human_v37_MT.fasta
```

Note that the index files need to remain in the same directory as the reference.
