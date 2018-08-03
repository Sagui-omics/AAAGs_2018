# Structural variants, BED format, and bedtools


## Structural variants (SVs)
Structural variants (SVs) are blocks of sequences that vary among individuals in their copy number (deletions, duplications, insertions), directionality (inversion), or location (translocation). The SVs have long been understudied, partly because of the difficulty of calling them. Currently there are four main methods for detecting SVs from NGS data, which are 1) read depth 2) split read 3) paired read, or 4) de novo assembly based. Many tools are developed using each of these methods, or combining multiple of them. However, there is still not a single perfect tool to call SVs with both fidelity and comprehensiveness. Therefore, when you receive a SV dataset, I would suggest you to understand the calling method and think about its limitation before you make a big statement out of a SV dataset.


## What is a BED format?

BED (Browser Extensible Data) format is commonly used to depict where the DNA fragments of interests are on a reference genome. In a BED file, each BED line represent a DNA fragment. The first three fields of a BED line tells you which chromosome this fragment is from, where it starts, and where it ends. The fields beyond these three are optional, you can keep whatever information it these fields, as long as all BED line has the same number of fields.

The BED coordinate is 0-based. For example, the first site on chr10 would be represent as: 

```
chr10	0	1
```

You can be downloaded the BED files of many different genomic features from [The Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)  of [UCSC Genome Browser](https://genome.ucsc.edu/).


## bedtools

bedtools is a powerful toolkit to analyze the spatial relationship between genomic fragments. It is commonly used in genome arithmetic such as overlapping, subtracting, and merging, of BED files. Note that some of the function are compatible with BAM and VCF file. With sufficient information, it can also convert a BED file to a BAM file, or vice versa. 

The [bedtools official site](http://bedtools.readthedocs.io/en/latest/#) has a tutorial. Or if you just want to quickly look up some of the functions, you can simply type in: 

```
bedtools --help
```
This will give you an overview of how the bedtools might be able to help you. 

Let us do it together. We start with activate the conda evironment so that we can use the bedtools installed using conda: 
```
export PATH=$PATH:/Users/sensnode/anaconda3/bin
source activate agar2018
```

Now again:

```
bedtools --help
```

## Demo 
### Determining and analyzing the hotspots of genomic structural variants 

The SVs do not occur uniformly across the human genome, but cluster into hotspots. I will lead you through the following steps to determine and to analyze the SV hotspots.

We will start with activate the conda evironment so that we can use the bedtools installed using conda.



### 1. Create interval windows across the human genome

First I am going to create a BED file with non-overlapping, 100Kb long intervals that cover the human reference genome so that I can count the number of SV(s) within them. This can be done using the `bedtools makewindow` function.
``` 
cd ~/Downloads/AAAGs_2018-master/Session7_StructuralVariants_Haplotypes/
bedtools makewindows -g hg19.genome -w 100000 > 100Kb_intervals.bed
```

### 2. Filter out the intervals overlapping with sequencing gaps

I want to exclude any intervals overlapping with a sequencing gap from the dataset because the SV number within these intervals are inevitably underestimated. This can be done using a subfunction of `bedtools intersect`. `bedtools intersect` basically report overlaps between two sets of intervals, specified after `-a` and `-b`. In this case, we use `-v` to make it report the intervals in the `-a` file that are NOT overlapping with any of the intervals in the `-b` file. 
```
bedtools intersect -a 100Kb_intervals.bed -b gap.bed -v > intervals_no_gap.bed
```

### 3. Count the number of SVs overlapping with each of the 100Kb intervals
To determine which of the intervals have unusually high number of overlaps with the SV(s), we want to first count how many SV(s) are there overlapping with each of the intervals. This can again be done using `bedtools intersect`. This time we use `-c`, which instruct bedtools to count the number of `-b` overlapping with each `-a`. The count number will be appended to the end of each `-a` line.   
```
bedtools intersect -a intervals_no_gap.bed -b SV_all_phase3.bed -c > intervals_SV.bed
```

Now I want to make a histogram out of the output file to see the distribution of the number of SV(s) each interval overlap with.

```
Rstudio
```
Just type `Rstudio` and it will open up the Rstudio and all the packages installed under our agar2018 environment.

In Rstudio, type in:
```
setwd("~/Downloads/AAAGs_2018-master/Session7_StructuralVariants_Haplotypes/")
t <- read.table("intervals_SV.bed", sep="\t")
names(t) <- c("chrom","pos","end","SVct")

library(ggplot2)
ggplot(t, aes(SVct)) + geom_histogram()
```

### 4. Create a mock SV dataset

Now I want to see what it would be like if the SV distribute randomly across the genome. I want to randomly permute the locations of the SVs . `bedtools shuffle` help create a mock set of intervals that are of the same number and sizes with the actual SVs set.   
```
bedtools shuffle -i SV_all_phase3.bed -g hg19.genome > mockSV.bed
```

### 5. Count the number of mock SVs overlapping with each of the 100Kb intervals
```
bedtools intersect -a intervals_SV.bed -b mockSV.bed -c > intervals_SV_mockSV.bed
```

We can make two overlaied histograms out of the file and see how the histogram of the SVs count in each interval is deviated from the histogram of the randomly distributed mock dataset.

```
t <- read.table("intervals_SV_mockSV.bed", sep="\t")
names(t) <- c("chrom","pos","end","SVct","mockSVct")

library(ggplot2)
ggplot(t) + 
	geom_histogram(aes(SVct),alpha=0.5)+
	geom_histogram(aes(mockSVct),fill="red",alpha=0.5)

```

### 6. Measure the coding sequence coverage of each of the 100Kb intervals
SVs are generally biased away from coding sequences. We therefore expect to see a lower level of coding sequences within the SV hotspots. Here we are going to measure the coding sequence coverage of each intervals. We use `bedtools coverage` here. This will report the depth and breadth of coverage of `-b` file on `-a` file. 
```
bedtools coverage -a intervals_SV_mockSV.bed -b CDS.bed > intervals_SV_mockSV_CDS.bed
```
With this output file, we can check if hotspots generally have lower level of coding sequences.

```
t <- read.table("intervals_SV_mockSV_CDS.bed", sep="\t")
names(t) <- c("chrom","pos","end","SVct","mockSVct", "CDSct", "CDSbp", "size","CDSfrac")

library(ggplot2)
ggplot(t, aes(SVct>5, CDSfrac, fill=SVct>5)) + geom_boxplot(notch=T, notchwidth=0.1, outlier.shape=NA)
```
