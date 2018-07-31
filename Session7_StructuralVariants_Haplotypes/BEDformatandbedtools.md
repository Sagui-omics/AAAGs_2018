# BED format and bedtools


## What is a BED format?

BED (Browser Extensible Data) format is commonly used to depict where the DNA fragments of interests are on a reference genome. In a BED file, each BED line represent a DNA fragment. The first three fields of a BED line tells you which chromosome this fragment is from, where it starts, and where it ends. The BED format is 0-based. For example, the first site on chr10 would be represent as: 

```
chr10	0	1
```



## bedtools

bedtools is a powerful toolkit to analyze the spatial relationship between genomic fragments. It is commonly used in genome arithmetic such as overlapping, subtracting, and merging, of BED files. Note that some of the function are compatible with BAM and VCF file. With sufficient information, it can also convert a BED file to a BAM file, or vice versa.

```
Yen-Lungs-iMac:~ onta$ bedtools -help
```
This will give you an overview of how the bedtools might be able to help you. I will show you some of them in the following session.



## Demo 
### Determining and analyzing the hotspots of genomic structural variants 

Structural variants (SVs) are blocks of sequences that vary among individuals in their copy number, directionality, or location. The SVs do not occur uniformly across the human genome, but cluster into hotspots.

I will lead you through the following steps to determine and to analyze the SV hotspots.



### 1. Create interval windows across the human genome

First I am going to create a BED file with non-overlapping, 100Kb long intervals that cover the human reference genome so that I can count the number of SV(s) within them. This can be done using the `bedtools makewindow` function.
``` 
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

### 4. Create a mock SV dataset

To see if the uneven distribution of SVs is random, I want to permute the locations of the SVs and compare the distribution. `bedtools shuffle` help create a mock set of intervals that are of the same number and sizes with the actual SVs set.   
```
bedtools shuffle -i SV_all_phase3.bed -g hg19.genome > mockSV.bed
```

### 5. Count the number of mock SVs overlapping with each of the 100Kb intervals
```
bedtools intersect -a intervals_SV.bed -b mockSV.bed -c > intervals_SV_mockSV.bed
```

We can make two simple histograms out of the file and see how the SV distribution is deviated from the random distribution.


### 6. Measure the coding sequence coverage of each of the 100Kb intervals
SVs are generally biased away from coding sequences. We therefore expect to see a lower level of coding sequences within the SV hotspots. Here we are going to measure the coding sequence coverage of each intervals. We use `bedtools coverage` here. This will report the depth and breadth of coverage of `-b` file on `-a` file. 
```
bedtools coverage -a intervals_SV_mockSV.bed -b CDS.bed > intervals_SV_mockSV_CDS.bed
```
With this output file, we can analyze if the intervals with more SVs (hotspots) generally have lower level of coding sequences.


