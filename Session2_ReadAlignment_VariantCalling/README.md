# AGAR 2018 Session 2: *The basics of read mapping and variant calling*

*Content by Tim Webster (Arizona State University and University of Utah), 2018. Available under a GNU General Public v3 license. Much of this content has been updated and adapted from a previous version of this activity (also by Webster) [that you can find here](https://github.com/thw17/ASU_BIO543_Genome_Assembly_2018).*

This session features a hands on introduction to the basics of read mapping, BAM file processing, and variant calling. We will also discuss various file types that you'll encounter along the way: FASTA, FASTQ, SAM/BAM, and VCF.

**Updates**
2018-08-02: Updated ``rule qualimap_per_sample`` in README and example_snakefile to reflect Snakemake 5.2's requirement that directory output targets are declared. Updated ``conda create`` command in "Setting Up" to force installation of Snakemake 5.2. Add exception to .gitignore for chimp_MT.fasta

**Table of Contents**

- Reference-based vs. *de novo* assembly
- Brief discussion of reproducibility and version control
- Setting up
- Reference genomes and the FASTA format
- Sequencing reads and the FASTQ format
- Building our pipeline
	- Preparing our reference genome
	- Inspecting read quality
	- Trimming reads
	- Mapping reads to the reference genome and processing alignments
	- Marking duplicate reads in SAM/BAM files
	- Exploring SAM/BAM files
	- Variant calling
- Running our pipeline
- Extending our pipeline
- Take home points

## Reference-based vs. *de novo* assembly

**Objectives: Understand the difference between reference-based and *de novo* assembly**

There are, in general, two main flavors of genome assembly. *De novo* assembly involves taking raw sequencing reads and piecing them together into a genome assembly using only the information contained in the reads and their metadata (e.g., the sequences themselves, insert sizes, etc.). While a number of de novo assemblers exist and there's a great deal of work being done to improve algorithms, lengthen sequencing reads, develop methods to increase insert sizes, etc., *de novo* assembly remains challenging, expensive (usually 100x or greater sequencing depth), and computationally demanding.

Fortunately, if we have a reference genome available to us, we can make do with much less sequencing (often 30x coverage or less; low coverage - 1-5x - is not uncommon for some purposes), and use tools that require far less memory and storage. We can do this by mapping our sequencing reads to said reference genome.

You can think of the two strategies as different ways of putting together a 30 million piece, highly repetitive jigsaw puzzle (assuming a 3 billion base pair genome and 100 base pair sequencing reads) with parts of many pieces clipped off (sequencing errors). In this analogy, *de novo* assembly is like mixing 100 or more of those puzzles together in a bag and then trying to put the puzzle together upside down. On the other hand, reference-guided assembly would be like mixing somewhere between 1 and 30 of these puzzles together and then putting the puzzle together picture up on top of a full-sized picture of what the puzzle is supposed to look like. Hopefully, from this analogy it should be clear that reference-guided assembly is much easier, provided your reference is close to what you're assembling and of good quality.

In this tutorial, we'll walk through the basics of reference-guided genome assembly. While the dataset we're working with is tiny (we're using the human mitochondrial genome and a tiny subset of reads from the 1000 genomes project), you should be able to use this as a starting point for working with larger datasets down the road.

**Questions for further thought/discussion**
1. When would we choose *de novo* assembly instead of reference-guided assembly?

2. In the case of reference-guided assembly, what are the pros and cons of higher and lower sequencing coverage?

## Brief discussion of reproducibility and version control

**Objectives: Understand motivation behind reproducibility and options for designing reproducible pipelines**

The fields of genomics and bioinformatics have typically embraced open source data and code, and reproducibility. While not all analyses are reproducible, these (and more) fields are continuing to move in this direction, with more and more journals requiring data and code to be shared.

With that in mind, it's important that any bioinformatic pipeline we put together can be shared in a way that allows someone to exactly (or as close as possible) reproduce our analyses. Because our time in this session is limited, a detailed discussion of reproducibility is beyond the scope of this tutorial. However, we will be building a reproducible, version controlled pipeline as we go using Anaconda/Bioconda (software manager) and Snakemake.

#### Conda and Bioconda

[Anaconda (conda)](https://conda.io/docs/) is an environment and package manager for the programming language Python and it makes installation, environment management, etc. simple without requiring root or administrator privileges. Fortunately, its framework has been leveraged to manage a variety of other languages and programs, including for a project called Bioconda that extends these capabilities to external bioinformatic programs as well. You can find out more in the [Bioconda documentation](https://bioconda.github.io/), the [Bioconda paper](https://www.nature.com/articles/s41592-018-0046-7), and the [Conda documentation](https://conda.io/docs/). We'll use Conda and Bioconda to create a controlled environment in which we can manage software specific to this pipeline. When we're done, we can use Conda to print the contents of this environment (every piece of software installed, plus its version and source; [details here](https://conda.io/docs/user-guide/tasks/manage-environments.html#sharing-an-environment)) so we can share it. In general, **I recommend creating a new Conda environment for every project.**

#### Snakemake

We'll be writing the pipeline itself using [Snakemake](https://snakemake.readthedocs.io/en/stable/), a workflow management system based on GNU Make, but written in Python and with added features that are specifically designed to aid in bioinformatic analyses. The [documentation](https://snakemake.readthedocs.io/en/stable/) is great, and they have [an excellent tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) that I highly recommend working through.

As I mentioned above, we won't be giving a Conda and/or Snakemake tutorial, but we will be using both. I explain our Conda set up below in the section **Setting up**. I'll very briefly describe our Snakemake setup here.

By default, Snakemake will read through your "Snakefile" (main file containing Snakemake code) until it hits the first "rule". This becomes the target rule for the entire pipeline. Snakemake will read the files required as input for this target rule, and then go through the output of the rest of the rules until it finds the rules that can create the input for this target rule. It then finds the rules whose output is the input of these rules, and so on until it has the complete path of rules (from the beginning) required to create the input for the target rule.

By convention, we'll call our first rule ``rule all`` and give the rest informative names. An example of a very simple Snakefile would be:

```
rule all:
	input:
		"hello.txt"

rule write_hello:
	output:
		"hello.txt"
	shell:
		"echo 'hello!!' > {output}"

rule write_goodbye:
	output:
		"goodbye.txt"
	shell:
		"echo 'goodbye!!' > {output}"
```
The first rule listed in the file is ``rule all``, so Snakemake checks it to see what input is required. Under ``input``, it sees that it needs to make ``hello.txt``. It then scans through the file, checking the output of each rule until it finds ``hello.txt``. As you can see, this is the case in ``rule write_hello``. ``rule write_hello`` has no input, so Snakemake will then run the command listed under ``shell`` in this rule (``echo 'hello!! > hello.txt``) and quit. Note that this command uses the command-line tool ``echo`` to write the string ``hello!!`` to the file ``hello.txt``.

There are three things to note in this example. 1) The command written in the ``shell`` section of the rule is the exact command that will be run on the command line (i.e., a UNIX style command). 2) Notice the ``{output}`` part of the shell command. This takes whatever is written in the ``output`` section of the rule and inserts it into the command. So, in ``rule write_hello``, the string under ``echo 'hello!!' > {output}`` becomes ``echo 'hello!!' > hello.txt``. 3) ``rule write_goodbye`` is completely ignored because its output, ``goodbye.txt``, isn't needed by ``rule all`` or any other rule required to generate the input of ``rule all`` (in this case, that's just ``rule write_hello``).

What if we want to write *both* ``hello.txt`` and ``goodbye.txt``? One option is to make the input of ``rule all`` a list, and just add ``goodbye.txt``. This would look like:

```
rule all:
	input:
		"hello.txt",
		"goodbye.txt"

rule write_hello:
	output:
		"hello.txt"
	shell:
		"echo 'hello!!' > {output}"

rule write_goodbye:
	output:
		"goodbye.txt"
	shell:
		"echo 'goodbye!!' > {output}"
```

**Note the comma in the input of ``rule all``**

This option works well in our case. Like before, Snakemake will run through the file until it hits the first rule (``rule all``). It now sees that it needs to make two files. For the first, ``hello.txt``, like before, it will see that this file is output by ``rule write_hello`` and knows that it has to run this rule. It then sees that ``goodbye.txt`` is output by ``rule write_goodbye`` and will run this rule. So, both rules will be run.

But what if we have many files that we want to output? We could list them all, but that can get quite tedious and might lead to human error. To help make this process easier (you'll see it become very useful in our pipeline below), Snakemake has a function called ``expand()`` that takes a string with brackets indicating text that needs to be expanded, and then the values to expand the brackets with. This should hopefully become clear in the following example. Take a look at ``rule all`` now.

```
rule all:
	input:
		expand("{word}.txt", word=["hello", "goodbye"])

rule write_hello:
	output:
		"hello.txt"
	shell:
		"echo 'hello!!' > {output}"

rule write_goodbye:
	output:
		"goodbye.txt"
	shell:
		"echo 'goodbye!!' > {output}"
```

In ``rule all``, Snakemake sees ``{word}.txt`` in the ``expand()`` function in ``input``. It then knows to replace the bracketed text with every item in the list called ``word`` after the comma (but within the parentheses of ``expand``). This results in a two item list for ``input``: ``hello.txt`` and ``goodbye.txt``.

Finally, to run a Snakemake pipeline, you can type (assuming Snakemake has been installed correctly and added to your PATH: see below):

``snakemake --snakefile <name of your snakefile>``

If your "Snakefile" is named ``Snakefile`` (note the capitalization), you can run your pipeline by just typing:

``snakemake``

Finally, to see what Snakemake is planning to run, you can run the command:

``snakemake -np``

which will just print what Snakemake is planning to do, without actually running the commands (i.e., a dry-run).

That's it for our Snakemake crash course for today. Obviously, we haven't begun to scratch the surface of Snakemake's capabilities or built an understanding of how it works. For this, I highly, highly recommend checking out the Snakemake [tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) and [documentation](https://snakemake.readthedocs.io/en/stable/).

**Questions for further though/discussion**
1. What are benefits of reproducibility for you and the community?

2. Why is it important to set up a new Conda environment for each project?

## Setting up

**Objectives: Clone git repository for this tutorial and use Conda to set up software environment for our pipeline**

For today's tutorial, you'll need this repository and Conda.

#### A quick note on notation

In this tutorial, all commands entered will follow a ```$ ```, while any associated output will follow this line in the
same box.  Do not include the ```$``` in your command, rather enter the command that follows.

#### Getting the repo
We'll use ```git``` to clone the repository (repo) for this tutorial onto your computer.  You can check to see if you have git installed by typing the command ```$ git```.  You should see some usage information.  If not, see [here](https://git-scm.com/) for information about installing.

Once git is installed, move to the directory on your computer or cluster where you'd like to work and type the command:

```
$ git clone https://github.com/Sagui-omics/AAAGs_2018
```
This should create a directory called AAAGs_2018 containing all of the workshop's materials.

Alternatively, if git isn't working for you, you can directly download the zipped directory via the green "Clone or download" button on [the repository's website](https://github.com/Sagui-omics/AAAGs_2018).

Now we need to change into our directory for this tutorial (all subsequent commands assume you're starting from here):

```
$ cd AAAGs_2018/Session2_ReadAlignment_VariantCalling
```

#### Setting up Anaconda
We're going to use [Conda](https://conda.io/docs/), described above, to install and manage software. To download Conda and set up our environment, use the following steps:

* First, install Python 3.6 version of Miniconda [available here](https://conda.io/miniconda.html), *OR* Anaconda [available here](https://www.continuum.io/downloads). Miniconda installs the conda framework without all of the Python packages installed with Anaconda (numpy, scipy, matplotlib, etc.). All of these packages can be easily installed later (via ``conda install <package name>``), so the decision is up to you. During installation, be sure to allow Miniconda/Anaconda to append to your .bashrc or .bash_profile (this will add it and all programs it installs to your PATH). *This means you'll have to pay attention to all prompts during installation!!*  If installation goes well, the command ``` which python ``` should result in something like ``/Users/<yourusername>/miniconda/bin/python `` or ``/home/<yourusername>/miniconda/bin/python ``.

* Add Bioconda channels to conda with the following commands in this order:

  ```
  $ conda config --add channels defaults

  $ conda config --add channels conda-forge

  $ conda config --add channels bioconda
  ```

	This sets the channel ``bioconda`` as our highest priority, followed by ``conda-forge``, and then ``defaults`` at the lowest priority. You can confirm this information (and more) with the command:

  ```
  $ conda info
  ```

	And you can manually add, delete, and adjust channels in the file ``~/.condarc``. This is a hidden file (you can tell because it begins with a single period) located in your home directory (denoted by the ``~/``). However, I don't recommend editing this file until you're familiar with both conda and UNIX/LINUX environments.

* Create the environment we'll be working in and install required packages with the command:

  ```
  $ conda create -n agar2018 python=3.6 'snakemake>=5.2 samtools bwa bioawk fastqc multiqc bbmap qualimap gatk4 vcftools picard
  ```

	This will create a working environment called agar2018 containing python 3.6 (python 3 is required for snakemake) and all of the tools listed in the command.  You can see the full list of programs available through Bioconda [listed here](https://bioconda.github.io/) and the full list of python packages available through Anaconda [listed here](https://docs.continuum.io/anaconda/pkg-docs). Note that there are many other channels to check as well (particularly, conda-forge).

If you want to load our new environment, you can do so with the command:
```
$ source activate agar2018
```
and leave the environment with the command:
```
$ source deactivate
```

If you're in your environment, you can easily add additional programs (like we did for samblaster) and packages with the command:
```
$ conda install <program/package name>
```

For example, if we also want to take a look at ``Bowtie2``, another read mapper (we'll use ``bwa`` today), we can easily add it by entering our environment ``` $ source activate agar2018 ``` and entering the command ```$ conda install bowtie2 ```

## Reference genomes and the FASTA format

**Objectives: Understand the FASTA format**

Because we're using a reference-guided assembly approach in this tutorial, we need a reference genome to which we're going to map reads. In a perfect world, this assembly is of a high-quality, has a good set of annotations available (e.g., genes, functional elements, repeats, etc.), and is relatively closely related to the species that you're studying.  There are, for example, numerous reference genomes hosted at the [UCSC Genome Browser](http://hgdownload.soe.ucsc.edu/downloads.html) and [Ensembl](http://www.ensembl.org/info/data/ftp/index.html).  Accessing a reference genome probably won't be a problem if you're working with a model organism, but in other situations you'll have to consider whether a good assembly is available for a taxon evolutionarily close enough for your purposes (if not, you might need to think about assembling a reference for your project _de novo_).  For today, we're working with example sequencing reads from human samples and we have the human reference genome available, so we'll be fine.

FASTA format is a file format designed for DNA or protein sequences.  It looks something like (the first 10 lines of the ``` human_v37_MT.fasta``` file in the ```references``` directory:

```
$ head reference/human_g1k_v37_MT.fasta

>MT
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
CGTCTGGGGGGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTC
GCAGTATCTGTCTTTGATTCCTGCCTCATCCTATTATTTATCGCACCTACGTTCAATATT
ACAGGCGAACATACTTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA
ACAATTGAATGTCTGCACAGCCACTTTCCACACAGACATCATAACAAAAAATTTCCACCA
AACCCCCCCTCCCCCGCTTCTGGCCACAGCACTTAAACACATCTCTGCCAAACCCCAAAA
ACAAAGAACCCTAACACCAGCCTAACCAGATTTCAAATTTTATCTTTTGGCGGTATGCAC
TTTTAACAGTCACCCCCCAACTAACACATTATTTTCCCCTCCCACTCCCATACTACTAAT
CTCATCAATACAACCCCCGCCCATCCTACCCAGCACACACACACCGCTGCTAACCCCATA
```
FASTA format includes a single ID line starting with a ``>``, followed by sequence lines. This repeats for each sequence in the file (one ID line, any number of sequence lines) with a space in between each sequence. In our example, the name of the sequence, ``MT`` is given after ``>`` on the ID line.  The lines that follow contain the sequence of ``MT``.  Because most FASTA files wrap lines every 50-80 characters (this isn't uniform across files, unfortunately), there will often be many lines composing each sequence.  Today's (```human_v37_MT.fasta```) file should only contain a single sequence, the 1000 genomes reference MT sequence, that's a bit more than 16 kb (i.e., 16 thousand base pairs) long.  We can quickly check to make sure using (the very, very powerful) [bioawk](https://github.com/lh3/bioawk), which we installed earlier:

```
$ bioawk -c fastx '{print ($name), length($seq)}' reference/human_v37_MT.fasta

MT	16569
```

We see that we do indeed have a single sequence called "MT" that's 16,569 bases in length.

## Sequencing reads and the FASTQ format

**Objectives: Understand the FASTQ format**

There are a few types of sequencing technologies out there, but Illumina is still dominant, so that's what we'll focus on today.  We won't get into the nuts and bolts of sequencing itself (there are plenty of resources available online that describe it well like [this video](https://www.youtube.com/watch?v=fCd6B5HRaZ8) or [this review](http://www.nature.com/nrg/journal/v17/n6/abs/nrg.2016.49.html)).

If you've sent your samples to a core for sequencing, they'll likely return to you a series of FASTQ files.  If you used single-end sequencing, your files can be concatenated into a single FASTQ file per sample per lane (note that you can easily concatenate files, even if they are gzipped, using ```cat```). On the other hand, if you used paired-end sequencing, you'll end up with two FASTQ files per sample per lane - one for the forward read and one for the reverse read (see the video I linked to in the previous paragraph for more information about paired-end sequencing).  It's generally very important that paired reads are in the same order in both files.  If you're getting reads directly from a sequencing center, they should already organized this way.  However, you might have to sort and re-pair reads if you have, for example, stripped reads from a BAM file containing an alignment (the README in the ```fastq``` directory explains how to do this, if needed).

As an example, I have provided paired-end reads from one male sifaka in the ``example_files`` directory.  Forward reads are in ``example_files/example_read1.fastq.gz`` and reverse reads are in ``example_files/example_read2.fastq.gz``. The ``.gz`` extension indicates that I have compressed the files. We can take a look at the first two reads in the file ``example_files/example_read1.fastq.gz`` using ``zcat`` on Linux or ``gzcat`` on Mac/Unix combined with ``head``. ``zcat`` and ``gzcat`` will decompress the file contents on the fly.
  ```
$ zcat example_files/example_read1.fastq.gz | head -n 8

@HWI-ST0831:187:D1W6GACXX:4:1101:1428:1980 1:N:0:TTACCATGACCA
NTTTCGGGTTACATCCCAACAAATACTACAAACACTCACAAGGCAAGATGTTTACATATCGATTTTTTTCCTTGTT
+
#4=DFFFDFDFHAHIJIGIJHIIIIJJHIICHGHGIIJGIFGGGH>GHGH=@CFAFGIGGEDEHIGHFDCCEDDC@
@HWI-ST0831:187:D1W6GACXX:4:1101:1583:1912 1:N:0:TGACCATGACCA
NGGGTGGACCCAGTGTTTCTGTTAGTGGAAGAAGCGGTTCTGAAGAAAATGATTGAATTTATTGGCTGGGAAGAAG
+
#1=DADDFHHHHHCFGIJJJJJJIJIJIHHHGHIIHI@DHIJHHIIDGGGHHHFGJFAHIJGHIJEHHH#######
```
In these files, each sequencing read is listed in a series of four lines:
* the sequence identifier line, beginning with @
* the sequence line (consisting of A, T, C, G, or N)
* a comment line (here, the comment lines only contain +)
* a quality score line (ASCII characters)

The sequence identifier can contain a lot of information (see [Illumina's description for more information](http://support.illumina.com/content/dam/illumina-support/help/BaseSpaceHelp_v2/Content/Vault/Informatics/Sequencing_Analysis/BS/swSEQ_mBS_FASTQFiles.htm)), the combination of which will identify individual reads uniquely.  In this example, IDs are organized as: instrumentID:instrument_run:flowcellID:lane:tile:x_position:y_position.

The sequence line here is 76 bases long and contains the nucleotide sequence corresponding to each read.

The comment line usually either is a lone + or a + and then the sequence identifier repeated.  You'll generally ignore this line, but it's good to be aware of what's on it in case you're using basic shell command-line tools to do things like count.  For example, if the sequence identifier is repeated and you're counting the number of times a tile number is present in a file using a command like ```grep```, your count will be double the actual number of reads coming from that tile.

The quality score line contains an ASCII character for every nucleotide in the sequence (i.e., it'll be 250 characters long for a 250 base read, 75 characters long for a 75 base read, etc.).  Each ASCII character can be converted into a integer PHRED quality score, ranging from 0 to 93, that indicates the probability that a particular base call is incorrect.  [See more details here - our Illumina scores use the 33 offset](http://www.drive5.com/usearch/manual/quality_score.html).

Like FASTA files, we can also use ```bioawk``` to parse and analyze fastq files.  For example we can count all of the reads in each file:

```
$ for i in example_files/*.fastq.gz; do echo $i; bioawk -c fastx 'END{print NR}' $i; done

example_files/example_read1.fastq.gz
4
example_files/example_read2.fastq.gz
4
```

As you can see, there are four reads in each file. This is an example of a "for loop" in BASH.  The ```echo``` command is simply printing each file name and then the ```bioawk``` command counts and prints the total number of records in each file.

We can also count the number of reads from tile 1101 in ``example_files/example_read1.fastq.gz``:

```
$ bioawk -c fastx '{print $name}' example_files/example_read1.fastq.gz | grep ':1101:' | wc -l

4
```

Or count the number of reads from each tile in ```example_files/example_read1.fastq.gz```:

```
$ bioawk -c fastx '{print $name}' example_files/example_read1.fastq.gz | cut -d':' -f 5 | sort | uniq -c

4 1101

```

In both cases, we see that all four reads come from the same tile. This is because I grabbed the first four reads from much larger files, and coming out of the sequencer, reads are sorted by tile.

**Questions for further thought/discussion**
1. Why are quality scores encoded as ASCII characters?

2. Bioinformatic programs generally require paired-end reads to be two files (forward read file and reverse read file), and to have reads in the two files sorted in the same order. Why is this organization important?

3. Can you tell what sequencing strategy and data type are included in a raw FASTQ file (e.g. DNA, RNA, whole genome, exome, RADseq, etc.)?

## Building our pipeline

Now let's move on to building our pipeline. I'm assuming that everything is downloaded, set up, and loaded correctly (see "Setting up" above).

As I mentioned in the "Brief discussion of reproducibility and version control" section above, this pipeline will be using Snakemake and Conda, but we won't be discussing them in detail. Please consult the links I provided above (documentation, tutorials, etc.) for more information.

So, open up your favorite text editor (I prefer Atom, but others like Sublime, Text Wrangler/BBEdit, Vim/Vi, emacs, and nano; just don't use things like Microsoft Word). And let's get started!

Let's first open up a new file and save it as "Snakefile" in the main project directory.

Next, let's add some Python variables that will store our sample IDs and a shorthand name for our reference. Copy and paste the following into Snakefile and save it:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]
```

As you can hopefully see, today we're working with data from 18 samples, six each from three populations: CEU, PUR, and YRI.

So far we don't have any Snakemake code, but that'll soon change.

#### Preparing our reference genome

**Objectives: Index a reference genome using BWA and SAMtools**

Remember from our discussion of reference genomes above, that the FASTA format consists of, per sequence, an ID line and a series of sequence lines. Most reference genomes are quite large, so it's very inefficient to linearly search through, say, 3 billion characters spread across 25 (or MANY more) sequences.  So, many programs use various hashing/indexing strategies to more efficiently handle the data. We'll create two of the most commonly required index files - .dict and .fai - which will summarize the reference sequences.  We'll also create the required index for ``bwa``, our read mapper, which will allow ``bwa`` to quickly search for matches in the reference.  This is all we'll need for our purposes today, but check any additional tools you incorporate into your work down the line to see if they require additional indexing or processing of reference assemblies.  For example, each read mapper will likely require its own unique index, so your ```bwa``` indices won't work for, say, ```bowtie2``` or ``hisat2``.

Today our reference is the mitochondrial genome sequence from the 1000 Genomes Project reference genome. See the README in the ``reference`` directory for more information about how it was extracted.

From the main project directory, run the following three commands:

```
$ samtools faidx reference/human_v37_MT.fasta

$ samtools dict -o reference/human_v37_MT.dict reference/human_v37_MT.fasta

$ bwa index reference/human_v37_MT.fasta
```

This will be quick and require very little memory on our small reference, but the bwa indexing in particular can take much longer on a full, human-sized reference and require ~4-6 GB of memory.

And that's it.  All of our reference files and indices are now contained in our reference directory. Let's take a quick look at our ``.fai`` index (created by the first command with ``samtools faidx``) with ``cat``:

```
$ cat reference/human_v37_MT.fasta.fai
MT	16569	4	60	61
```

The columns in this output are sequence name, sequence length, offset (in bytes), number of bases per line, and number of bytes per line. Generally speaking, as users, we will only be interested in the first two columns. And, checking those values, it's reassuring that our reference FASTA contains only one sequence, called "MT", that's 16,569 bases long. If you remember from before, that's the exact value we got using ``bioawk``.

Now it's time to add this information to our Snakefile. Let's add variables with our tool paths (in case we want to update later), a rule all (to run the whole pipeline), and a rule that will run those three commands to prepare our reference. Add this information, so our Snakefile now looks like:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]

# Tool paths
bwa_path = "bwa"
samtools_path = "samtools"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies)

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")
```
Note that the variable ``assemblies`` in the input of ``rule all`` is the same ``assemblies`` declared above. It's a list of one item, ``"human_v37_MT"``, a string.

Save that file. From the main directory, what happens when we type the command to do a dry run of the pipeline? It should look like this:

```
$ snakemake -np
Building DAG of jobs...
Nothing to be done.
```

Why is there nothing to be done? Because we already created ``reference/human_v37_MT.fasta.fai`` manually when we tested our commands above. Let's force Snakemake to ignore what we've already done. We can do that with:

```
$ snakemake -np --forceall

Building DAG of jobs...
Job counts:
	count	jobs
	1	all
	1	prepare_reference
	2

rule prepare_reference:
    input: reference/human_v37_MT.fasta
    output: reference/human_v37_MT.fasta.fai, reference/human_v37_MT.fasta.amb, reference/human_v37_MT.dict
    jobid: 1
    wildcards: assembly=human_v37_MT


localrule all:
    input: reference/human_v37_MT.fasta.fai
    jobid: 0

Job counts:
	count	jobs
	1	all
	1	prepare_reference
	2
```

You can see that Snakemake is now going to run ``rule prepare_reference`` to create the input of ``rule all``.

**Questions for further thought/discussion**
1. Why did I include ``{assembly}`` as a variable, when we only have one value in ``assemblies``?

2. Take a look at the contents of ``reference/human_v37_MT.dict`` with the command ``cat reference/human_v37_MT.dict``. What information does it contain? What does M5 mean and why is that value important?

3. Generally speaking, you should map sequencing reads to the whole genome and not just a single region (we're doing the latter today to allow analyses to be run quickly on a local computer). Why do you think this is the case?

#### Inspecting read quality

**Objectives: Use fastqc and multiqc to explore sequencing read quality in FASTQ files**

The first thing you should do when you get FASTQ files is get a sense of their quality.  The quickest way to do this is to run ```fastqc``` and take a look at the reports.

We've provided 36 FASTQ files: forward and reverse reads from 18 individuals from three populations in the 1000 Genomes Project. These are located in the ``fastq`` directory. The README in the ``fastq`` directory has more information about the samples and how we obtained the files.

Let's generate a fastqc report for one of the samples. For this we can use the command:

```
$ fastqc -o fastqc_results fastq/CEU_NA07000_MT.*
```
The asterisk means "anything", so this command can be read as "use fastqc on anything that starts with fastq/CEU_NA07000_MT. and output to the directory called fastqc_results."

This command should take a few seconds to complete and output a .zip and .html file for each of the two FASTQs (it'll take longer for full-sized files).  You can open the .html locally on your computer in your browser. Let's open up both of the .html files.

Overall, they seem to be of good quality. Per base sequence quality stays high across the read (i.e., above 30). GC and sequence content look good. We also don't see much in the way of sequence duplication or overrepresented sequences.

If your sequences are of an unexpectedly low quality it might be worth contacting your sequencing center.  Otherwise, lower quality towards the ends of the reads, some PCR duplication, slightly lower quality on the reverse read, and sequence content that's slightly off are all pretty typical in sequencing experiments.  Also note that certain sequencing experiments, for example those targeting coding regions like exome sequencing and RNA seq, are likely to have some odd kmers.

What does a very bad report look like? Take a look at the example I provided: ``example_files/SRR740818_2_fastq.html``. This comes from a large public dataset and is a very good example of a failed run.

Obviously, it'd be quite tedious to run that command for each of our 18 samples. It's also hard to view the reports individually. So, let's use the power of Snakemake to automate our ``fastqc`` analyses and use a tool called ``multiqc`` to summarize the findings. If we add these rules, our file should now look like (be sure to save it):

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]

# Tool paths
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
multiqc_path = "multiqc"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html"

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_MT.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_MT.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_MT.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"
```
In addition to the two new rules, note the additions to the input of ``rule all``, and the new program path variables.

We can then run our full pipeline with the command:

```
$ snakemake
```

This should take a couple minutes to run (it has to generate 36 reports). When it's done, take a look at ``multiqc_results/multiqc_report.html``. You can see that it contains an interactive summary of all of the fastqc reports.

**Questions for furtherthought/discussion**
1. In the updated Snakefile, why didn't we include any of the output of ``rule fastqc_analysis`` as input to ``rule all``?

2. Why are we concerned about FASTQ quality? How might it affect downstream analyses?

#### Trimming reads

**Objectives: Remove adapters and trim low-quality ends of sequencing reads**

After getting an initial assessment of FASTQ quality. It's time to think about trimming your reads. Opinions on whether or not to trim vary immensely. In general, you should remove sequencing adapters if possible and then carefully consider whether or not to trim for quality. One thing to consider is that DNA mappers like BWA tend to handle adapters and low-quality sequence well, while many RNA mappers have more trouble. If you're at all concerned, I recommend trying a variety of trimming parameters to see how they affect your data.

To give you a sense of tools and commands that you might use, we're going to trim using ``bbduk.sh``. I've provided Illumina adapter sequences in ``misc/adapter_sequence.fa``. The following command will remove these adapter sequences as well as do a little bit of conservative quality trimming for sample ``PUR_HG00553``:

```
$ bbduk.sh -Xmx1g in1=fastq/PUR_HG00553_MT.R1.fastq.gz in2=fastq/PUR_HG00553_MT.R2.fastq.gz out1=trimmed_fastqs/PUR_HG00553_trimmed_read1.fastq.gz out2=trimmed_fastqs/PUR_HG00553_trimmed_read2.fastq.gz ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe qtrim=rl trimq=15 minlen=50 maq=20
```

We should then run ``fastqc`` and ``multiqc`` on the trimmed files so that we can compare the trimmed reports to our original reports.

If we add rules for these steps to our Snakefile, it'll now read:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]

# Tool paths
bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
multiqc_path = "multiqc"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html",
		"multiqc_trimmed_results/multiqc_report.html"

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_MT.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_MT.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_MT.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = "fastq/{sample}_MT.R1.fastq.gz",
		fq2 = "fastq/{sample}_MT.R2.fastq.gz"
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx1g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_results/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed_results/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed_results/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples, read=["read1", "read2"])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed_results fastqc_trimmed_results"
```

Open up both multiqc reports in a browser (``multiqc_trimmed_results/multiqc_report.html`` and ``multiqc_results/multiqc_report.html``) and see how trimming affected our data.

**Questions for further thought/discussion**
1. Did trimming affect any of the samples? Focus in particular on mean quality scores, N content, and adapter content.

2. Are there differences among population in sequence quality and subsequent trimming results?

**Alternative Programs**
While we focus on ``bbduk.sh`` here, I also recommend checking out [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html) and [Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/).

#### Mapping reads to the reference genome and processing alignments

**Objectives: map reads to reference genome and convert alignments to BAM format; understand SAM/BAM format**

Our next step involves mapping our reads to our reference.  We'll use the [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) algorithm to do this, as it's among the most popular mappers and works very well mapping reads to a closely related reference.  Other popular mappers include [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [minimap2](https://github.com/lh3/minimap2), [Novoalign](http://www.novocraft.com/products/novoalign/), and [Stampy](http://www.well.ox.ac.uk/stampy) (Stampy, in particular, for mapping to a very evolutionary diverged reference genome). Note that for RNA seq data you would use different mappers like [Star](https://github.com/alexdobin/STAR) and [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml), among many others.

The command line for ``bwa mem`` is quite straightforward.  From the main directory, we can execute the following command:

```
$ bwa mem reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz > bams/YRI_NA18498.sam
```

This will map both sets of paired-end reads from ``YRI_NA18498`` to our reference genome and output the alignments in SAM format in our ``bams`` directory.

If we take a look at the top of the file using ``$ head -n 4 bams/YRI_NA18498.sam``, we see the following:

```
@SQ	SN:MT	LN:16569
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz
SRR027524.6392421	83	MT	7053	60	76M	=	6904	-225	GGGGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCT	<<8?@B?@@CBCAADACCAAB@@BAABCBCBCCEABDAACDCAAAAAADCB>ECABCBBEABACBCAABE>>A>A>	NM:i:1	MD:Z:2A73	MC:Z:76M	AS:i:73	XS:i:0
SRR027524.6392421	163	MT	6904	60	76M	=	7053	225	GATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCCTGACTGGCATTGTATTAGC	;CA@BDABD@BB?D@C@CDBA?==C?B?>A>@BA@C>>@C;9;?D@38?@B6=A>A:B@B@BD@ACA=DA?>4<><	NM:i:0	MD:Z:76	MC:Z:76M	AS:i:76	XS:i:0
```

As you can see, our first line (@SQ) contains information about the reference genome (there would be more lines if there were more sequences in our reference), our second line (@PG) contains information about our ``bwa`` command, and the remaining lines contain information about each mapped read.  You can find more information about what information is contained in each record in the [SAM/BAM specificaions](https://samtools.github.io/hts-specs/SAMv1.pdf).

While it's exciting that we've successfully mapped our first sample, there are two potential problems with our command above.  First, while SAM format is convenient in that it's human-readable, alignment files are ENORMOUS, so we'll want to compress (BAM format) to minimize our data footprint.  Further, we want to ensure that all read-pairing, etc. didn't get lost in the mapping process.  We can easily add these steps to our pipeline by piping our output to ```samtools``` which handles this easily.  So our complete command now looks like this:
```
$ bwa mem reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz | samtools fixmate -O bam - bams/YRI_NA18498.bam
```

Great! We've mapped reads, corrected read pairing, and converted to BAM format! However, the bam file we just created contains all of our alignments, but it's currently unordered, unlabeled, unfiltered, and unindexed.

##### Adding read groups

**Objectives: add read group information to a BAM file, learn how to view a BAM file**

Read groups are very useful when you're working with multiple samples, sequencing lanes, flowcells, etc.  Importantly for us, it'll help our downstream variant caller label samples correctly and handle potential sequencing batch effects.  [Picard](https://broadinstitute.github.io/picard/command-line-overview.html) is very commonly used to add read groups to bam files, but ```bwa``` also has the ability to add read groups on the fly while mapping.  This latter option will save us time and space, so we'll add read groups to individual ``YRI_NA18498`` with ```bwa``` by adding to our previous command:
  ```
  $ bwa mem -M -R '@RG\tID:YRI_NA18498\tSM:YRI_NA18498\tLB:YRI_NA18498\tPU:YRI_NA18498\tPL:Illumina' reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz | samtools fixmate -O bam - bams/YRI_NA18498.bam
  ```
In a perfect world, we'd know more about our sample and could use these tags more appropriately. For now we're including them as placeholders that you can fill in for your future pipelines. ID is the name of a read group (containing a unique combination of the following tags).  SM is the sample name.  LB is the sequencing library. PU is the flowcell barcode.  PL is the sequencing technology.  A number of other options exist (see the [@RG section of the SAM/BAM specifications](https://samtools.github.io/hts-specs/SAMv1.pdf) for more information).

*Note that the FASTQ files we're using today actually contain reads from multiple lanes, etc., as we randomly grabbed them from high-coverage 1000 genomes bam files.  But for simplicity's sake in this tutorial, we'll ignore that.*

Let's take a look at our new BAM file. Because BAM files are compressed, we can't view them using our standard command line tools. Instead, we'll use ``samtools view``. Here's a command to print the first 5 lines of the file:

```
$ samtools view -h bams/YRI_NA18498.bam | head -n 5

@SQ	SN:MT	LN:16569
@RG	ID:YRI_NA18498	SM:YRI_NA18498	LB:YRI_NA18498	PU:YRI_NA18498	PL:Illumina
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -M -R @RG\tID:YRI_NA18498\tSM:YRI_NA18498\tLB:YRI_NA18498\tPU:YRI_NA18498\tPL:Illumina reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz
SRR027524.6392421	83	MT	7053	60	76M	=	6904	-225	GGGGCTGTATTTGCCATCATAGGAGGCTTCATTCACTGATTTCCCCTATTCTCAGGCTACACCCTAGACCAAACCT	<<8?@B?@@CBCAADACCAAB@@BAABCBCBCCEABDAACDCAAAAAADCB>ECABCBBEABACBCAABE>>A>A>	NM:i:1	MD:Z:2A73	AS:i:73	XS:i:0	RG:Z:YRI_NA18498	MQ:i:60MC:Z:76M
SRR027524.6392421	163	MT	6904	60	76M	=	7053	225	GATCTGCTGCAGTGCTCTGAGCCCTAGGATTCATCTTTCTTTTCACCGTAGGTGGCCTGACTGGCATTGTATTAGC	;CA@BDABD@BB?D@C@CDBA?==C?B?>A>@BA@C>>@C;9;?D@38?@B6=A>A:B@B@BD@ACA=DA?>4<><	NM:i:0	MD:Z:76	AS:i:76	XS:i:0	RG:Z:YRI_NA18498	MQ:i:60	MC:Z:76M
```

Notice the new ``@RG`` line (2nd line) that gives the read group ID (``ID``) and the associated values. You should also be able to see that each read now has a ``RG`` tag (``RG:Z:YRI_NA18498``) telling us what read group that particular read belongs to. In our example we only have one read group, but you can imagine how this becomes useful if a BAM file contains data from multiple sequencing lanes, sequencing libraries, samples, etc.

##### Sorting bam files

**Objectives: learn how to sort a BAM file in genomic coordinate order**

Sorting doesn't require too much explanation.  Most genomic datasets are huge, so it's inefficient to move along unsorted bam files (we need access to all reads covering a given base for variant calling, for example).  ```samtools sort ``` is widely used, and that's what we'll employ here.  Like our previous tools, it handles streaming input, so we can simply add to our previous command to save space:
  ```
  $ bwa mem -M -R '@RG\tID:YRI_NA18498\tSM:YRI_NA18498\tLB:YRI_NA18498\tPU:YRI_NA18498\tPL:Illumina' reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz | samtools fixmate -O bam - - | samtools sort -O bam -o bams/YRI_NA18498.sorted.bam -
  ```

Let's take a quick look at our new, sorted BAM using ``samtools view`` again.

```
$ samtools view -h bams/YRI_NA18498.sorted.bam | head -n 5

@HD	VN:1.5	SO:coordinate
@SQ	SN:MT	LN:16569
@RG	ID:YRI_NA18498	SM:YRI_NA18498	LB:YRI_NA18498	PU:YRI_NA18498	PL:Illumina
@PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:bwa mem -M -R @RG\tID:YRI_NA18498\tSM:YRI_NA18498\tLB:YRI_NA18498\tPU:YRI_NA18498\tPL:Illumina reference/human_v37_MT.fasta trimmed_fastqs/YRI_NA18498_trimmed_read1.fastq.gz trimmed_fastqs/YRI_NA18498_trimmed_read2.fastq.gz
SRR027523.7329778	353	MT	1	60	44H32M	=	83	154	GATCACAGGTCTATCACCCTATTAACCACTCA	AB>@DADB@@ACAA@D@@@BAA;AB@:B>@=NM:i:0	MD:Z:32	MC:Z:72M	AS:i:32	XS:i:0	RG:Z:YRI_NA18498	SA:Z:MT,16526,+,44M32S,60,0;
```

You can see that there is one major change to the header.  There is a @HD line in the header indicating that the file is sorted in coordinate order ("SO:coordinate"),

##### Indexing

Again, BAM files can get pretty big and they're compressed, so we need to index them for other tools to use them. Here's a simple command for indexing our bam (from our main directory):
```
$ samtools index bams/YRI_NA18498.sorted.bam
```

This creates a file with a .bai extension that needs to remain in the same directory as its corresponding bam. Also, we're going to have to index every BAM file we create as we move through the pipeline.

Now we're at a good point to update our Snakefile and run our mapping and BAM processing steps up to this point. After adding our new rules for mapping and indexing, and updating ``rule all``, our Snakefile will now read:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]

# Tool paths
bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
multiqc_path = "multiqc"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html",
		"multiqc_trimmed_results/multiqc_report.html",
		expand(
			"bams/{sample}.{assembly}.sorted.bam.bai",
			sample=all_samples,
			assembly=assemblies)

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_MT.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_MT.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_MT.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = "fastq/{sample}_MT.R1.fastq.gz",
		fq2 = "fastq/{sample}_MT.R2.fastq.gz"
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx1g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_results/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed_results/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed_results/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples, read=["read1", "read2"])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed_results fastqc_trimmed_results"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai"
	output:
		"bams/{sample}.{assembly}.sorted.bam"
	params:
		id = "{sample}",
		sm = "{sample}",
		lb = "{sample}",
		pu = "{sample}",
		pl = "Illumina",
		bwa = bwa_path,
		samtools = samtools_path
	shell:
		" {params.bwa} mem -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_sorted_bam:
	input:
		"bams/{sample}.{assembly}.sorted.bam"
	output:
		"bams/{sample}.{assembly}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"
```

Save this command and run it using:

```
$ snakemake
```

We can see what files were output in the ``bams`` directory with the command:

```
$ ls bams/
CEU_NA06984.human_v37_MT.sorted.bam
CEU_NA06984.human_v37_MT.sorted.bam.bai
CEU_NA06985.human_v37_MT.sorted.bam
CEU_NA06985.human_v37_MT.sorted.bam.bai
CEU_NA06986.human_v37_MT.sorted.bam
CEU_NA06986.human_v37_MT.sorted.bam.bai
CEU_NA06989.human_v37_MT.sorted.bam
CEU_NA06989.human_v37_MT.sorted.bam.bai
CEU_NA06994.human_v37_MT.sorted.bam
CEU_NA06994.human_v37_MT.sorted.bam.bai
CEU_NA07000.human_v37_MT.sorted.bam
CEU_NA07000.human_v37_MT.sorted.bam.bai
PUR_HG00551.human_v37_MT.sorted.bam
PUR_HG00551.human_v37_MT.sorted.bam.bai
PUR_HG00553.human_v37_MT.sorted.bam
PUR_HG00553.human_v37_MT.sorted.bam.bai
PUR_HG00554.human_v37_MT.sorted.bam
PUR_HG00554.human_v37_MT.sorted.bam.bai
PUR_HG00637.human_v37_MT.sorted.bam
PUR_HG00637.human_v37_MT.sorted.bam.bai
PUR_HG00638.human_v37_MT.sorted.bam
PUR_HG00638.human_v37_MT.sorted.bam.bai
PUR_HG00640.human_v37_MT.sorted.bam
PUR_HG00640.human_v37_MT.sorted.bam.bai
YRI_NA18486.human_v37_MT.sorted.bam
YRI_NA18486.human_v37_MT.sorted.bam.bai
YRI_NA18488.human_v37_MT.sorted.bam
YRI_NA18488.human_v37_MT.sorted.bam.bai
YRI_NA18489.human_v37_MT.sorted.bam
YRI_NA18489.human_v37_MT.sorted.bam.bai
YRI_NA18498.human_v37_MT.sorted.bam
YRI_NA18498.human_v37_MT.sorted.bam.bai
YRI_NA18499.human_v37_MT.sorted.bam
YRI_NA18499.human_v37_MT.sorted.bam.bai
YRI_NA18501.human_v37_MT.sorted.bam
YRI_NA18501.human_v37_MT.sorted.bam.bai
```

A sorted BAM and corresponding index are present for every sample.

#### Marking duplicate reads in BAM files

**Objectives: introduce duplicate reads and how to mark/remove them**

During library preparation for sequencing, amplification steps can lead to PCR duplicates of reads.  The inclusion of duplicate reads can negatively affect our downstream analyses, so we need to mark or remove them (marking allows downstream tools to ignore them). Much like trimming, opinions differ on how important it is to remove duplicates. With improving library preparation methods, duplicates probably only have a minor effect on downstream variant calling. However, in cases of lower quality samples, of targeted capture, or when count data matters, duplication will have a greater effect. As long as you're working with smaller datasets and you're not computationally limited, it safe enough to mark/remove duplicates.

Today, we're going to use [Picard](https://broadinstitute.github.io/picard/command-line-overview.html), as it's going to be required if you have a single sample sequenced across multiple sequencing flowcells. If you're working with smaller datasets where each sample is only sequenced on a single flowcell, I recommend taking a look at  [Samblaster](https://github.com/GregoryFaust/samblaster), which can take streaming output from ``bwa``, speeding up duplicate removal and allowing us to produce one less bam file (saving us space). [Sambamba](http://lomereiter.github.io/sambamba/) is another nice option.  Note that at the time of writing this tutorial, [the developers of Samtools suggest against using it to remove duplicates](https://github.com/samtools/samtools/issues/614).

Using ``picard`` we can mark duplicates with the following command:

```
$ picard -Xmx1g MarkDuplicates I=bams/PUR_HG00640.human_v37_MT.sorted.bam O=bams/PUR_HG00640.human_v37_MT.sorted.mkdup.bam M=stats/PUR_HG00640.human_v37_MT.picard_mkdup_metrics.txt
```
This will flag, but not remove, duplicates in our bam.

A quick note here.  If the same sample library is sequenced across multiple lanes, you'll want to identify and mark/remove duplicates AFTER merging the bam files from all of the lanes. We're not merging today, however.

Remember that we have to index these BAM files as well. After adding these two steps (marking duplicates and indexing) to the Snakefile, it will look like:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]

# Tool paths
bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
multiqc_path = "multiqc"
picard_path = "picard"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html",
		"multiqc_trimmed_results/multiqc_report.html",
		expand(
			"bams/{sample}.{assembly}.sorted.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"bams/{sample}.{assembly}.sorted.mkdup.bam.bai",
			sample=all_samples,
			assembly=assemblies),

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_MT.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_MT.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_MT.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = "fastq/{sample}_MT.R1.fastq.gz",
		fq2 = "fastq/{sample}_MT.R2.fastq.gz"
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx1g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_results/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed_results/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed_results/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples, read=["read1", "read2"])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed_results fastqc_trimmed_results"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai"
	output:
		"bams/{sample}.{assembly}.sorted.bam"
	params:
		id = "{sample}",
		sm = "{sample}",
		lb = "{sample}",
		pu = "{sample}",
		pl = "Illumina",
		bwa = bwa_path,
		samtools = samtools_path
	shell:
		" {params.bwa} mem -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_sorted_bam:
	input:
		"bams/{sample}.{assembly}.sorted.bam"
	output:
		"bams/{sample}.{assembly}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "bams/{sample}.{assembly}.sorted.bam",
		bai = "bams/{sample}.{assembly}.sorted.bam.bai"
	output:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		metrics = "stats/{sample}.{assembly}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx1g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"bams/{sample}.{assembly}.sorted.mkdup.bam"
	output:
		"bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"
```

Run the additional steps of the pipeline with:

```
$ snakemake

```

**Questions for further thought/discussion**
1. Why do you think Picard can only remove duplicates from a BAM file in genomic coordinate sorted order?

2. If duplicates can possibly affect even one variant, why is there debate about removing them? What are possible arguments against duplicate removal?

#### Exploring SAM/BAM files

**Objectives: calculate various BAM file metrics**

We've already seen how we can view BAM files using ``samtools view``, but there are usually millions of sequencing reads. How can we efficiently summarize a BAM file?

[Samtools offers tools to calculate summary statistics](http://www.htslib.org/doc/samtools.html).  Today, we'll calculate stats using ```samtools stats``` because it provides a bit more detail, but you should have a look at ```samtools flagstat``` as well.  From our main directory, enter the command:

```
$ samtools stats bams/CEU_NA06994.human_v37_MT.sorted.mkdup.bam | grep ^SN | cut -f 2- > stats/CEU_NA06994.human_v37_MT.sorted.mkdup.bam.stats
```

Because samtools stats offers a huge range of statistics including a number of very big tables in our output, we'll just grab the summary statistics.  This is what ```grep ^SN | cut -f 2-``` in our command does.

We can print the contents of each file to screen using the ```cat``` command.  Here's the command and its result for the samtools stats output for ```stats/CEU_NA06994.human_v37_MT.sorted.mkdup.bam.stats```:

```
cat stats/CEU_NA06994.human_v37_MT.sorted.mkdup.bam.stats
raw total sequences:	9376
filtered sequences:	0
sequences:	9376
is sorted:	1
1st fragments:	4688
last fragments:	4688
reads mapped:	9322
reads mapped and paired:	9268	# paired-end technology bit set + both mates mapped
reads unmapped:	54
reads properly paired:	8962	# proper-pair bit set
reads paired:	9376	# paired-end technology bit set
reads duplicated:	77	# PCR or optical duplicate bit set
reads MQ0:	0	# mapped and MQ=0
reads QC failed:	0
non-primary alignments:	0
total length:	919410	# ignores clipping
bases mapped:	914075	# ignores clipping
bases mapped (cigar):	911968	# more accurate
bases trimmed:	0
bases duplicated:	7559
mismatches:	2556	# from NM fields
error rate:	2.802730e-03	# mismatches / bases mapped (cigar)
average length:	98
maximum length:	100
average quality:	32.9
insert size average:	388.4
insert size standard deviation:	866.3
inward oriented pairs:	4582
outward oriented pairs:	59
pairs with other orientation:	2
pairs on different chromosomes:	0
```

As you can see, it gives us a lot of information about number of reads, mapping, pairing, etc.  A quick glance shows us that our mapping was quite successful (9322 out of 9376 reads mapped).  We also had very few PCR duplicates (77 reads - note that this can be calculated because we flagged, but didn't remove duplicates), but this is probably because we randomly sampled read pairs to create our FASTQ files.

We can also get some additional statistics and helpful summary plots and figures using ``qualimap``:

```
$ qualimap bamqc -bam bams/CEU_NA06986.human_v37_MT.sorted.mkdup.bam -nt 1 -outdir stats/qualimap/human_v37_MT/CEU_NA06986/
```

This command writes a number of files, but we're interested in the one ending in ``.html``, which we can open up in a browser.

We can now add these steps (along with a step for Qualimap to summarize from all samples together) to the Snakefile. After doing so, our Snakefile will look like:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]

# Tool paths
bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
multiqc_path = "multiqc"
picard_path = "picard"
qualimap_path = "qualimap"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html",
		"multiqc_trimmed_results/multiqc_report.html",
		expand(
			"bams/{sample}.{assembly}.sorted.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"bams/{sample}.{assembly}.sorted.mkdup.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"stats/{sample}.{assembly}.sorted.mkdup.bam.stats",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"stats/qualimap/{assembly}/multisampleBamQcReport.html",
			assembly=assemblies)

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_MT.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_MT.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_MT.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = "fastq/{sample}_MT.R1.fastq.gz",
		fq2 = "fastq/{sample}_MT.R2.fastq.gz"
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx1g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_results/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed_results/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed_results/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples, read=["read1", "read2"])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed_results fastqc_trimmed_results"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai"
	output:
		"bams/{sample}.{assembly}.sorted.bam"
	params:
		id = "{sample}",
		sm = "{sample}",
		lb = "{sample}",
		pu = "{sample}",
		pl = "Illumina",
		bwa = bwa_path,
		samtools = samtools_path
	shell:
		" {params.bwa} mem -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_sorted_bam:
	input:
		"bams/{sample}.{assembly}.sorted.bam"
	output:
		"bams/{sample}.{assembly}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "bams/{sample}.{assembly}.sorted.bam",
		bai = "bams/{sample}.{assembly}.sorted.bam.bai"
	output:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		metrics = "stats/{sample}.{assembly}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx1g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"bams/{sample}.{assembly}.sorted.mkdup.bam"
	output:
		"bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		"stats/{sample}.{assembly}.sorted.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule qualimap_per_sample:
	input:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		directory("stats/qualimap/{assembly}/{sample}")
	params:
		qualimap = qualimap_path,
		out_dir = "stats/qualimap/{assembly}/{sample}/"
	shell:
		"{params.qualimap} bamqc "
		"-bam {input.bam} -nt 1 "
		"-outdir {params.out_dir}"


rule create_qualimap_list:
	input:
		lambda wildcards: expand(
			"stats/qualimap/{genome}/{sample}",
			sample=all_samples,
			genome=wildcards.assembly)
	output:
		"stats/qualimap/{assembly}/qualimap.list"
	run:
		shell("echo -n > {output}")
		for i in input:
			sm = os.path.basename(i)
			shell("echo '{}\t{}' >> {{output}}".format(sm, i))

rule qualimap_multibamqc:
	input:
		"stats/qualimap/{assembly}/qualimap.list"
	output:
		"stats/qualimap/{assembly}/multisampleBamQcReport.html"
	params:
		qualimap = qualimap_path,
		out_dir = "stats/qualimap/{assembly}/"
	shell:
		"{params.qualimap} multi-bamqc -d {input} -outdir {params.out_dir}"

```

This will take a couple minutes to run. When it's done, take a look at the Qualimap file summarizing all of the samples, ``stats/qualimap/human_v37_MT/multisampleBamQcReport.html``, which you can open in your browser.

**Questions for further thought/discussion**
1. We started with 10,000 reads per sample (5,000 forward and 5,000 reverse)--why are there fewer reads in each of the BAM files?

2. Is coverage uniform across the mitochondrial genome? How does this affect our target depth of coverage for sequencing?

#### Variant calling

**Objectives: jointly call variants across all samples**

Now that we've trimmed and mapped our reads, labeled and sorted reads in our BAM files, removed duplicates, and assessed our mapping success, it's finally time to call variants. By "call variants", I mean statistically infer genotypes.

There are a few tools that do this (see Alternative Programs at the end of this section), but today we're going to use the [Genome Analysis Toolkit's (GATK's) HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.7.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php). This involves a three-step process of first preliminarily genotyping each sample individually, then combining the all of these preliminary genotypes into a single file, and finally jointly genotyping across all samples. This has a few benefits. 1) It scales well across many, many samples. 2) Joint genotyping can increase the power to identify difficult alleles. For example, if a certain allele is present in a small fraction of reads for individual 1, it might be ignored if individual 1 is called by itself. However, if that allele is present in individual 3, it can increase the support for that allele in individual 1. This can be extremely important for low-coverage sequencing. 3) GATK's HaplotypeCaller locally reassembles regions of the genome. This isn't unique to HaplotypeCaller and it's benefits are beyond the scope of this tutorial, but this helps with things like small insertions and deletions that can cause reads to locally mismap.

For step one of this process, we use HaplotypeCaller to preliminarily genotype each sample separately:

```
$ gatk --java-options "-Xmx1g" HaplotypeCaller -R reference/human_v37_MT.fasta -I bams/CEU_NA07000.human_v37_MT.sorted.mkdup.bam -ERC GVCF -O gvcfs/CEU_NA07000.human_v37_MT.g.vcf.gz
```
Note the option ``-ERC GVCF``, which outputs confidence that an invariant site (homozygous and matches reference genome) is a reference allele. This is important for the joint genotyping.

Our next step is run once the previous step is finished for all samples. It uses GATK's CombineGVCFs to combine the GVCF files from all samples. You'll see what this command looks like in the Snakefile in a minute.

The final step involves jointly genotyping that combined GVCF file. For this we use GATK's GenotypeGVCFs tool, and again, you'll see what the command looks like in the Snakefile.

Let's add these final rules to our Snakefile. Our final Snakefile, containing our full pipeline from start to finish, looks like:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT"]

# Tool paths
bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
gatk_path = "gatk"
multiqc_path = "multiqc"
picard_path = "picard"
qualimap_path = "qualimap"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html",
		"multiqc_trimmed_results/multiqc_report.html",
		expand(
			"bams/{sample}.{assembly}.sorted.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"bams/{sample}.{assembly}.sorted.mkdup.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"stats/{sample}.{assembly}.sorted.mkdup.bam.stats",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"stats/qualimap/{assembly}/multisampleBamQcReport.html",
			assembly=assemblies),
		expand(
			"genotyped_vcfs/{assembly}.gatk.called.raw.vcf.gz",
			assembly=assemblies)

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_MT.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_MT.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_MT.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = "fastq/{sample}_MT.R1.fastq.gz",
		fq2 = "fastq/{sample}_MT.R2.fastq.gz"
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx1g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_results/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed_results/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed_results/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples, read=["read1", "read2"])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed_results fastqc_trimmed_results"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai"
	output:
		"bams/{sample}.{assembly}.sorted.bam"
	params:
		id = "{sample}",
		sm = "{sample}",
		lb = "{sample}",
		pu = "{sample}",
		pl = "Illumina",
		bwa = bwa_path,
		samtools = samtools_path
	shell:
		" {params.bwa} mem -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_sorted_bam:
	input:
		"bams/{sample}.{assembly}.sorted.bam"
	output:
		"bams/{sample}.{assembly}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "bams/{sample}.{assembly}.sorted.bam",
		bai = "bams/{sample}.{assembly}.sorted.bam.bai"
	output:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		metrics = "stats/{sample}.{assembly}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx1g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"bams/{sample}.{assembly}.sorted.mkdup.bam"
	output:
		"bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		"stats/{sample}.{assembly}.sorted.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule qualimap_per_sample:
	input:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		directory("stats/qualimap/{assembly}/{sample}")
	params:
		qualimap = qualimap_path,
		out_dir = "stats/qualimap/{assembly}/{sample}/"
	shell:
		"{params.qualimap} bamqc "
		"-bam {input.bam} -nt 1 "
		"-outdir {params.out_dir}"


rule create_qualimap_list:
	input:
		lambda wildcards: expand(
			"stats/qualimap/{genome}/{sample}",
			sample=all_samples,
			genome=wildcards.assembly)
	output:
		"stats/qualimap/{assembly}/qualimap.list"
	run:
		shell("echo -n > {output}")
		for i in input:
			sm = os.path.basename(i)
			shell("echo '{}\t{}' >> {{output}}".format(sm, i))

rule qualimap_multibamqc:
	input:
		"stats/qualimap/{assembly}/qualimap.list"
	output:
		"stats/qualimap/{assembly}/multisampleBamQcReport.html"
	params:
		qualimap = qualimap_path,
		out_dir = "stats/qualimap/{assembly}/"
	shell:
		"{params.qualimap} multi-bamqc -d {input} -outdir {params.out_dir}"

rule gatk_gvcf:
	input:
		ref = "reference/{assembly}.fasta",
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		"gvcfs/{sample}.{assembly}.g.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} """
		"""-ERC GVCF -O {output}"""

rule gatk_combinegvcfs:
	input:
		ref = "reference/{assembly}.fasta",
		gvcfs = lambda wildcards: expand(
			"gvcfs/{sample}.{genome}.g.vcf.gz",
			sample=all_samples,
			genome=wildcards.assembly)
	output:
		"combined_gvcfs/{assembly}.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		print(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}"""
		)
		shell(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf:
	input:
		ref = "reference/{assembly}.fasta",
		gvcf = "combined_gvcfs/{assembly}.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs/{assembly}.gatk.called.raw.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""

```

You can run the remaining steps of the pipeline by again typing:

```
$ snakemake
```

The final output file, ``genotyped_vcfs/human_v37_MT.gatk.called.raw.vcf.gz`` is an unfiltered VCF file containing only variant sites. Viewing, understanding, and filtering this file is the subject of the next session, taught by Maria Nieves Colon.

**Questions for further thought/discussion**
1. Are there times that it might be inappropriate to jointly genotype samples?

**Additional Programs**
I recommend also checking out [Freebayes](https://github.com/ekg/freebayes). If you're working with low coverage data, you should consider working with [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

## Running our pipeline
As we went along today, I hope you saw that Snakemake didn't rerun previously completed commands. This is because Snakemake checks for for input/output files. It won't rerun a rule if the output is present and newer than the input. If, however, you happened to update the timestamp of the reference genome:

```
$ touch reference/human_v37_MT.fasta
```

and print a dry run of the pipeline:

```
$ snakemake -np
```

You should see that Snakemake plans to rerun all of the rules except for FASTQ quality assessment and read trimming. This is because these rules do not depend on the reference genome.

[Snakemake has many, many options worth reading about in the documentation](https://snakemake.readthedocs.io/en/stable/), but a few handy commands for running your pipeline are:

1. Run your pipeline that's written in a file called "Snakemake"
```
$ snakemake
```

2. Run your pipeline that's written in a file called "something_else"
```
$ snakemake --snakefile something_else
```

3. Print the steps that Snakemake plans to take, but don't run anything
```
$ snakemake -np
```

4. Print the steps that Snakemake plans to take, don't run anything, and tell us why each rule needs to be run
```
$ snakemake -npr
```

5. Run your full pipeline from the beginning, no matter what's been done already

```
$ snakemake --forceall
```

Again, check the Snakemake documentation and tutorial for more information. Snakemake also has excellent support for deploying a pipeline across a computing cluster and all of that information is available in the documentation.

## Extending our pipeline

Once you get the hang of Snakemake, you'll realize that it's quite easy to extend your pipelines.

For example, throughout the pipeline, I've been including a variable for assembly in file names, but we've been using a list of assemblies with only one value. This is by design because I work on both genome assembly and with nonmodel taxa, so I frequently map to multiple genome versions. Let's quickly add a new genome to the pipeline.

Note that ``reference`` also contains a chimpanzee mitochondrial reference genome: ``reference/chimp_MT.fasta``. To add this genome to the pipeline and run things in parallel, mapping to both the human and chimp references, simply add ``"chimp_MT"`` to the assemblies list. Our full Snakefile now reads:

```
ceu = [
	"CEU_NA06984",
	"CEU_NA06985",
	"CEU_NA06986",
	"CEU_NA06989",
	"CEU_NA06994",
	"CEU_NA07000"]

pur = [
	"PUR_HG00551",
	"PUR_HG00553",
	"PUR_HG00554",
	"PUR_HG00637",
	"PUR_HG00638",
	"PUR_HG00640"]

yri = [
	"YRI_NA18486",
	"YRI_NA18488",
	"YRI_NA18489",
	"YRI_NA18498",
	"YRI_NA18499",
	"YRI_NA18501"]

all_samples = ceu + pur + yri
assemblies = ["human_v37_MT", "chimp_MT"]

# Tool paths
bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
samtools_path = "samtools"
fastqc_path = "fastqc"
gatk_path = "gatk"
multiqc_path = "multiqc"
picard_path = "picard"
qualimap_path = "qualimap"

rule all:
	input:
		expand(
			"reference/{assembly}.fasta.fai",
			assembly=assemblies),
		"multiqc_results/multiqc_report.html",
		"multiqc_trimmed_results/multiqc_report.html",
		expand(
			"bams/{sample}.{assembly}.sorted.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"bams/{sample}.{assembly}.sorted.mkdup.bam.bai",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"stats/{sample}.{assembly}.sorted.mkdup.bam.stats",
			sample=all_samples,
			assembly=assemblies),
		expand(
			"stats/qualimap/{assembly}/multisampleBamQcReport.html",
			assembly=assemblies),
		expand(
			"genotyped_vcfs/{assembly}.gatk.called.rawvariants.vcf.gz",
			assembly=assemblies)

rule prepare_reference:
	input:
		ref = "reference/{assembly}.fasta"
	output:
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell(
			"{params.samtools} faidx {input.ref}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {input.ref}")
		# bwa
		shell(
			"{params.bwa} index {input.ref}")

rule fastqc_analysis:
	input:
		"fastq/{sample}_MT.{read}.fastq.gz"
	output:
		"fastqc_results/{sample}_MT.{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_results {input}"

rule multiqc_analysis:
	input:
		expand(
			"fastqc_results/{sample}_MT.{read}_fastqc.html",
			sample=all_samples,
			read=["R1", "R2"])
	output:
		"multiqc_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f "
		"-o multiqc_results fastqc_results"

rule trim_adapters_paired_bbduk:
	input:
		fq1 = "fastq/{sample}_MT.R1.fastq.gz",
		fq2 = "fastq/{sample}_MT.R2.fastq.gz"
	output:
		out_fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx1g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=50 maq=20"

rule fastqc_analysis_trimmed:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		html1 = "fastqc_trimmed_results/{sample}_trimmed_read1_fastqc.html",
		html2 = "fastqc_trimmed_results/{sample}_trimmed_read2_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_results {input.fq1} {input.fq2}"

rule multiqc_analysis_trimmed:
	input:
		expand(
			"fastqc_trimmed_results/{sample}_trimmed_{read}_fastqc.html",
			sample=all_samples, read=["read1", "read2"])
	output:
		"multiqc_trimmed_results/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} --interactive -f -o multiqc_trimmed_results fastqc_trimmed_results"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_fastqs/{sample}_trimmed_read2.fastq.gz",
		ref = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai"
	output:
		"bams/{sample}.{assembly}.sorted.bam"
	params:
		id = "{sample}",
		sm = "{sample}",
		lb = "{sample}",
		pu = "{sample}",
		pl = "Illumina",
		bwa = bwa_path,
		samtools = samtools_path
	shell:
		" {params.bwa} mem -R "
	 	"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_sorted_bam:
	input:
		"bams/{sample}.{assembly}.sorted.bam"
	output:
		"bams/{sample}.{assembly}.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule picard_mkdups:
	input:
		bam = "bams/{sample}.{assembly}.sorted.bam",
		bai = "bams/{sample}.{assembly}.sorted.bam.bai"
	output:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		metrics = "stats/{sample}.{assembly}.picard_mkdup_metrics.txt"
	params:
		picard = picard_path
	shell:
		"{params.picard} -Xmx1g MarkDuplicates I={input.bam} O={output.bam} "
		"M={output.metrics}"

rule index_mkdup_bam:
	input:
		"bams/{sample}.{assembly}.sorted.mkdup.bam"
	output:
		"bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats:
	input:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		"stats/{sample}.{assembly}.sorted.mkdup.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule qualimap_per_sample:
	input:
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		directory("stats/qualimap/{assembly}/{sample}")
	params:
		qualimap = qualimap_path,
		out_dir = "stats/qualimap/{assembly}/{sample}/"
	shell:
		"{params.qualimap} bamqc "
		"-bam {input.bam} -nt 1 "
		"-outdir {params.out_dir}"


rule create_qualimap_list:
	input:
		lambda wildcards: expand(
			"stats/qualimap/{genome}/{sample}",
			sample=all_samples,
			genome=wildcards.assembly)
	output:
		"stats/qualimap/{assembly}/qualimap.list"
	run:
		shell("echo -n > {output}")
		for i in input:
			sm = os.path.basename(i)
			shell("echo '{}\t{}' >> {{output}}".format(sm, i))

rule qualimap_multibamqc:
	input:
		"stats/qualimap/{assembly}/qualimap.list"
	output:
		"stats/qualimap/{assembly}/multisampleBamQcReport.html"
	params:
		qualimap = qualimap_path,
		out_dir = "stats/qualimap/{assembly}/"
	shell:
		"{params.qualimap} multi-bamqc -d {input} -outdir {params.out_dir}"

rule gatk_gvcf:
	input:
		ref = "reference/{assembly}.fasta",
		bam = "bams/{sample}.{assembly}.sorted.mkdup.bam",
		bai = "bams/{sample}.{assembly}.sorted.mkdup.bam.bai"
	output:
		"gvcfs/{sample}.{assembly}.g.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} """
		"""-ERC GVCF -O {output}"""

rule gatk_combinegvcfs:
	input:
		ref = "reference/{assembly}.fasta",
		gvcfs = lambda wildcards: expand(
			"gvcfs/{sample}.{genome}.g.vcf.gz",
			sample=all_samples,
			genome=wildcards.assembly)
	output:
		"combined_gvcfs/{assembly}.gatk.combinegvcf.g.vcf.gz"
	params:
		gatk = gatk_path
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		print(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}"""
		)
		shell(
			"""{params.gatk} --java-options "-Xmx1g" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf:
	input:
		ref = "reference/{assembly}.fasta",
		gvcf = "combined_gvcfs/{assembly}.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs/{assembly}.gatk.called.rawvariants.vcf.gz"
	params:
		gatk = gatk_path
	shell:
		"""{params.gatk} --java-options "-Xmx1g" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""
```

Note the *only* difference is the ``assemblies`` variable declaration near the top of the file. If we run ``snakemake -np``, we can see that Snakemake in fact does plan to do everything in parallel.

This is true for adding samples, running multiple mappers (e.g., a ``bwa`` vs. ``minimap2`` comparision), etc. As long as you name things and declare variables carefully, these type of extensions are very easy.

In addition, I declare tool names at the top of the file in case I or someone else has to use a different version of a tool (e.g., a local version of ``samtools`` that you downloaded). In this example, we only have to change ``samtools_path`` once at the top of the file, rather than in every single rule that uses ``samtools``.
