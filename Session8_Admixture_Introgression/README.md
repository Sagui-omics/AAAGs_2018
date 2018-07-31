# AGAR - Admixture and Genetic Introgression Analysis  

In this session you will learn how to perform analysis to infer introgression and admixture using the [ADMIXTOOLS] (https://github.com/DReichLab/AdmixTools) package. For more details on the theory being D-statistics and F3 tests see [Patterson et al. Ancient Admixture in Human History (2012) Genetics. vol. 192 no. 3 1065-1093](http://www.genetics.org/content/192/3/1065)

- You should aim to spend the first hour learning how to convert from VCF format to the PLINK and EIGENSTRAT formats.
- You should aim to spend the second hour learning how to perform a D-statistic analysis to test for introgression.
- You should aim spend the third how learning how to identify evidence of admixture using the f3 test.

The following exercises should be performed on the LINUX cluster provided. You have been provided a ***AGAR\_admixture.tar.gz*** file for this session. When you unpack it, the following files should be present. If any are missing, please let the instructor know.

File Contents | 
----------------
AltaiNea.vcf.gz |
AltaiNea.vcf.gz.tbi |
Haak.Yamnaya.bed |
Haak.Yamnaya.bim |
Haak.Yamnaya.fam |
Laz.ASW_MXL.geno |
Laz.ASW_MXL.ind |
Laz.ASW_MXL.snp |
map_pops.py |
pop_list |

##Section 1. Converting to Eigenstrat/Ancestry format

You have been provided a VCF file for the [high coverage Altai Neanderthal sequence] (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4031459/) called at a specific set of ~600K SNPs (***AltaiNea.vcf.gz***). This is a specially edited subset of the original Neanderthal data, which can be found [here] (cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF). Along with the vcf is its index file, ending ***.tbi*** .
Our aim is to get this into a biallelic SNP format that is usable by ADMIXTOOLS. There are a few ways to do this, but in this session we will first convert the VCF to a PLINK format, and then to an EIGENSTRAT format. We will not be doing analysis on the actual PLINK files, but these are usually useful to have regardless, as PLINK provides much more [flexibility for editing and filtering SNP data](https://www.cog-genomics.org/plink2) than offered by ADMIXTOOLS.

We must first convert the VCF file to PLINK format. The following command using [vcftools] (http://vcftools.sourceforge.net/man_latest.html) will do this.

	$ vcftools --gzvcf AltaiNea.vcf.gz --remove-indels --plink --out AltaiNea

- **vcftools** invokes the program vcftools
- **--gzvf** tells vcftools our input VCF file
- **--remove-indels** tells vcftools to ignore any indels in the VCF file.
- **--plink**' tells vcftools to convert the genotypes in the VCF file into a PLINK file<sup>1</sup>
- **--out** tells vcftools what the output prefix should be


This command should produce three files:
 ***AltaiNea.log***, ***AltaiNea.ped***, ***AltaiNea.map***

The ***.ped*** and ***.map*** files are our PLINK files. The former contains information about the individual samples and the and individual genotypes, and the latter contains information about the SNP positions. For more information see [here] (http://zzz.bwh.harvard.edu/plink/data.shtml#ped). These two files together are considered *flat* PLINK files. However it is common practice to compress them into a much more easily readable (and smaller) *binary* format. We do this using the following command:

	$ plink2 --file AltaiNea --make-bed --out AltaiNea

- **plink2** invokes the plink command
- **--file** tells plink what the input file is (it automatically looks for the ***ped*** and ***map*** parts of the file name)
- **--make-bed** tells plink to make a binary file
- **--out** tells plink what the output prefix should be

The command should produce three files:
 ***AltaiNea.bed***, ***AltaiNea.bim***, ***AltaiNea.fam***

The ***.bed*** file contains the genotype information, and is not human readable. The ***.bim*** file contains the snp info and the ***.fam*** file contains the individual info. These latter files two are human readable, and we should be able to find the name of the sample of interest, *AltaiNea*, in the ***.fam*** file. 

PLINK files are useful for a lot of things. But ADMIXTOOLS generally wants a slightly different format. The PLINK binary file in the ADMIXTOOLS manual is referred to as PACKEDPED format. We want to convert this to PACKEDANCESTRYMAP<sup>2</sup>. To do this, we use a tool found in ADMIXTOOLS called [convertf](https://github.com/DReichLab/AdmixTools/tree/master/convertf).

Like all ADMIXTOOLS programs, convertf requires a *parameter* file to work. In this file we provide the input files (in this case the PLINK files), the output filetype we want, and the names of the output files. convertf automatically works out the format of the input files (partly by the naming). We need to create the parameter file in an editor, such as **vi**.

	$ vi par_convertf_AltaiNea

The text we enter into the parameter file should look like this:

	genotypename: AltaiNea.bed
	snpname: AltaiNea.bim
	indivname: AltaiNea.fam
	outputformat: PACKEDANCESTRYMAP
	genotypeoutname: AltaiNea.geno
	snpoutname: AltaiNea.snp
	indivoutname: AltaiNea.ind

Now we can run convertf like so:

	$ convertf -p par_convertf_AltaiNea

- **convertf** invokes the convertf command
- **-p par\_convertf\_AltaiNea** tells convertf to read the ***par\_convertf\_AltaiNea*** parameter file.

As specified in our parameter file, this will result in the following output files.
***AltaiNea.geno***, ***AltaiNea.ind***, ***AltaiNea.fam***

Though the format is slightly different, these three files essentially contain the same information as the **.bed**, **.fam** and **.bim** PLINK files. The one weakness with PLINK for the kind of analysis ADMIXTOOLS does is it does not maintain *population information* associated with each sample. The current **.ind** file will look like this:

	AltaiNea:AltaiNea U        ???
	
We need to remove the *???* and put in the population our sample belongs to (for example **Altai**). We could do this manually, or we could use **sed**.

	$ sed -i 's\???\Altai\g' AltaiNea.ind

- **sed** invokes the program sed
- **-i** tells sed to change the file in place
- **s\???\Altai\g** tells sed to change all instances of *???* in the file with *Altai*

	

We have now made a file that can be used by ADMIXTOOLS. The resulting **.ind** file will now look like this:

	AltaiNea:AltaiNea U        Altai

Of course, there is not much we  can do with this file on its own (it is just data for one genome). Therefore we are going to merge this file with another set of data in PACKEDANCESTRYMAP format. This dataset, with the prefix ***Laz.ASW\_MXL***, consists of a subset of data from [Lazaridis et al. 2014] (https://www.nature.com/articles/nature13673) plus African American (ASW) and Mexican Amerian (MLX) 1000 Genomes populations.

ADMIXTOOLS does not have the diverse options of PLINK for manipulating SNP, but it does allow fairly smart merging of two SNP datasets using the **mergit** program. Again we must define a parameter file, containing the names of the datasets we want to merge, and their outputs. Use the general structure below to create the appropriate parameter file to merge the ***AltaiNea*** dataset with the ***Laz.ASW.MXL*** dataset. We suggest using ***merge\_AltaiNed\_par*** as the name of the parameter file and ***Laz.ASW\_MXL.AltaiNea*** as the prefix for the output files.

	geno1: first input genotype file
	snp1:  first input snp file
	ind1:  first input indiv file 
	geno2: second input genotype file
	snp2:  second input snp file
	ind2:  second input indiv file 
	genooutfilename:   output genotype file
	snpoutfilename:    output snp file
	indoutfilename:    output indiv file

<!--- 
	$ vi merge_AltaiNed_par 

	geno1: Laz.ASW_MXL.geno
	snp1:  Laz.ASW_MXL.snp
	ind1:  Laz.ASW_MXL.ind
	geno2: AltaiNea.geno
	snp2:  AltaiNea.snp
	ind2:  AltaiNea.ind 
	genooutfilename:   Laz.ASW_MXL.AltaiNea.geno
	snpoutfilename:    Laz.ASW_MXL.AltaiNea.snp
	indoutfilename:    Laz.ASW_MXL.AltaiNea.ind 
--->

Once we have our proposed parameter file, we create the merged file like so. 
	
	$ mergeit -p  merge_AltaiNed_par

The program will run to completion<sup>3</sup> and we should now have some merged output files that we can do some analysis on. 



##Section 2. D-statistic Introgression Analysis

We will now learn how to calcultates D-statistics using the classic example of Neanderthals, Eurasians and sub-Saharan Africans.

Given the population ordering of (*pop1*, *pop2*, *pop3*, *pop4*), D-statistics allow the testing of the hypothesis of whether one of two sister tax (*pop1* and *pop2*) are more closely related to *pop3*. *pop4* is set as an outgroup to orientate whether a particular allele at biallelic SNP is the ancestral or derived variant. Noting **A** as ancestral and **B** as derived, the D-statistic test examines the ratio of the number of times *pop1* and *pop3* share the derived allele (**BABA**) versus *pop2* and *pop3* (**ABBA**). If there are more BABA than ABBA SNPs then the D-statistic will be **positive (+)**, while it will be **negative (-)** in the opposite case. If there is no difference between *pop1* and *pop2* versus *pop3*, then the expectation of the D-statistic is **0**. Significance of the D-statistics from 0 is assessed using a **Z-score**.

The program for conducting D-statistic tests is [qpDstat](https://github.com/DReichLab/AdmixTools/blob/master/README.Dstatistics). Once again, we must define a parameter file. However, in addition we must define a *population* file, where each line contains the names of 4 populations (space or tab seperated) and represents one particular D-statistic test (using the *pop1*, *pop2*, *pop3*, *pop4* order described above). The classic application of the D-statistic is testing whether Neaderthals are closer to non-Africans versus Africans. Within the ***Laz.ASW\_MXL*** file there is SNP data for the following modern populations:

+ sub-Saharan Africans
	+ Yoruba
	+  Mbuti
+ European 
	+ Sardinian 
	+ French
	+ Orcadian 
	+ Spanish
+ South Asian 
	+ Sindhi
	+ Burusho
+ East Asian 
	+ Han
	+ Japanese
+ Southeast Asian
	+ Bougainville
	+ Papuan
+ Middle Eastern
	+ Egyptian
	+ Palestinian
+ American
	+ Karitiana

In addition, we have provided data for **Chimp** to act as an outgroup. To apply the D-statistic to identfy the signal of Neanderthal introgression in non-Africans, we may make a file using vi called ***f4\_pop\_Altai*** containing the following 15 tests<sup>4</sup>:

	Sardinian	Yoruba	Altai	Chimp
	Sardinian	Mbuti	Altai	Chimp
	French	Yoruba	Altai	Chimp
	French	Mbuti	Altai	Chimp
	Han	Yoruba	Altai	Chimp
	Han	Mbuti	Altai	Chimp
	Japanese	Yoruba	Altai	Chimp
	Japanese	Mbuti	Altai	Chimp
	Sardinian	Han	Altai	Chimp
	French	Han	Altai	Chimp
	Sardinian	Japanese	Altai	Chimp
	French	Japanese	Altai	Chimp
	Sardinian	French	Altai	Chimp
	Han	Japanese	Altai	Chimp
	Yoruba	Mbuti	Altai	Chimp

We would then create a parameter file called ***f4\_par\_Altai*** containing the following information:

	genotypename: Laz.ASW_MXL.AltaiNea.geno 
	snpname: Laz.ASW_MXL.AltaiNea.snp
	indivname: Laz.ASW_MXL.AltaiNea.ind
	popfilename: f4_pop_Altai

We can then run the D-statistic tests like so:

	$ qpDstat -p f4_par_Altai > f4_res_Altai

When the tests are finished the results (or errors if something went wrong) will be in ***f4\_res\_Altai***. The first set of lines descibe details of the run.

	qpDstat: parameter file: f4_par_Altai
	### THE INPUT PARAMETERS
	##PARAMETER NAME: VALUE
	genotypename: Laz.ASW_MXL.AltaiNea.geno
	snpname: Laz.ASW_MXL.AltaiNea.snp
	indivname: Laz.ASW_MXL.AltaiNea.ind
	popfilename: f4_pop_Altai
	## qpDstat version: 701
	packed geno read OK
	end of inpack
	number of quadruples 15
	  0            Sardinian   20
	  1               French   20
	  2                  Han   20
	  3             Japanese   20
	  4               Yoruba   20
	  5                Mbuti   10
	  6                Altai    1
	  7                Chimp    1
	jackknife block size:     0.050
	snps: 594822  indivs: 112
	number of blocks for jackknife: 550
	nrows, ncols: 112 594822


The rows we are interested in begin with **result:**
	
	result:  Sardinian     Yoruba      Altai      Chimp      0.0249     7.463 18106  17228 594822 
	result:  Sardinian      Mbuti      Altai      Chimp      0.0199     5.440 18782  	18050 594822 
	result:     French     Yoruba      Altai      Chimp      0.0259     7.981 18154  17238 594822

There are then 9 columns of interest for each of these **result:** rows. The first four reiterate the test populations in the order we provided. The 5th provides the actual D-statistic. The 6th provides the Z score. The 7th, 8th and 9th columns are the number of BABA, ABBA and total number of SNPs considered respectively. The Z-score is generally the column we are most interested in. Typically scores above | 3 | (3 standard deviations removed from the mean) are considered significant in these kind of applications. Recall that the sign of the Z-score matters in terms of directionality of the greatest similarity of ***pop3***.

If everything has gone to plan, we should see that using Europeans or Asian populations as *pop1* and sub-Saharan Africans as *pop2* gives very high Z scores, suggesting that non-Africans are significantly closer to the Altai Neanderthals than sub-Saharan Africans. However, when setting non-Africans as both *pop1* and *pop2* (or Africans as *pop1* and *pop2*), the results are close to 0, suggesting that all non-Africans are pretty equidistant genetically to Neanderthals.

The other classic application of the D-statistic is showing that individuals from Southeast Asian are closer to Denisovans than other modern populations. We have included a high coverage **Denisovan** individual in the ***Laz.ASW\_MXL*** dataset. See if you can test this using the principles learnt above. Can you also show that Altai and Denisova are closer to each other than modern populations (i.e. that they form a clade)?

<!---  
	$ vi f4_par_Denisovan

	genotypename: Laz.ASW_MXL.AltaiNea.geno 
	snpname: Laz.ASW_MXL.AltaiNea.snp
	indivname: Laz.ASW_MXL.AltaiNea.ind
	popfilename: f4_pop_Denisovan


	$ vi f4_pop_Denisovan

	Sardinian	Papuan	Denisovan	Chimp
	Sardinian	Bougainville	Denisovan	Chimp
	French	Papuan	Denisovan	Chimp
	French	Bougainville	Denisovan	Chimp
	Han	Papuan	Denisovan	Chimp
	Han	Bougainville	Denisovan	Chimp
	Japanese	Papuan	Denisovan	Chimp
	Japanese	Bougainville	Denisovan	Chimp
	Yoruba	Papuan	Denisovan	Chimp
	Mbuti	Bougainville	Denisovan	Chimp
	Sardinian	Yoruba	Denisovan	Chimp
	Sardinian	Mbuti	Denisovan	Chimp
	French	Yoruba	Denisovan	Chimp
	French	Mbuti	Denisovan	Chimp
	Han	Yoruba	Denisovan	Chimp
	Han	Mbuti	Denisovan	Chimp
	Japanese	Yoruba	Denisovan	Chimp
	Japanese	Mbuti	Denisovan	Chimp
	Sardinian	Han	Denisovan	Chimp
	French	Han	Denisovan	Chimp
	Sardinian	Japanese	Denisovan	Chimp
	French	Japanese	Denisovan	Chimp
	Sardinian	French	Denisovan	Chimp
	Han	Japanese	Denisovan	Chimp
	Yoruba	Mbuti	Denisovan	Chimp

	$ qpDstat -p f4_par_Denisovan > f4_res_Denisovan
--->


##Section 3. F3 Admixture Analysis

The other common analysis that is done with ADMIXTOOLS is to use the F3 test to look for significant evidence of admixture. Theory shows that given allele frequencies at a given SNP for three populations, *A*, *B* and *C*, the F3 statistic [(C-A)(C-B)] is only negative<sup>5</sup> if a target population (*C*) was the product of admixture between two source populations *A* and *B* (or their ancestors). This can be examined in ADMIXTOOLS using the program [qp3Pop](https://github.com/DReichLab/AdmixTools/blob/master/README.3PopTest).

A good example of the application of the F3 test is showing that modern European populations are a mix of three ancestral populations. Two of these are Paleolithic/Mesolithic hunter-gatherers and Neolithic farmers. We have included one hunter-gatherer genome (**Loschbour**) and one farmer (**LBK**) in our ***Laz.ASW\_MXL*** dataset. The setup for these tests is very similar to D-statistic tests. Again we need a parameter file (we suggest to call it ***f3\_par\_PaleoNeo***) and population file containg the three-way tests (***f3\_pop\_PaleoNeo***). We will perform the test by comparing all modern populations to our two potential source populations. The parameter file should look like so:

<!---$ vi f3_par_PaleoNeo --->

	genotypename: Laz.ASW_MXL.AltaiNea.geno 
	snpname: Laz.ASW_MXL.AltaiNea.snp
	indivname: Laz.ASW_MXL.AltaiNea.ind
	popfilename: f3_pop_PaleoNeo

<!---$ vi f3_pop_PaleoNeo --->

The population file should look like this. Again, there should be one test per line.

	Loschbour	LBK	Yoruba
	Loschbour	LBK	Mbuti
	Loschbour	LBK	Sardinian
	Loschbour	LBK	French
	Loschbour	LBK	Orcadian
	Loschbour	LBK	Spanish
	Loschbour	LBK	Sindhi
	Loschbour	LBK	Egyptian
	Loschbour	LBK	Palestinian
	Loschbour	LBK	Karitiana
	Loschbour	LBK	Burusho
	Loschbour	LBK	Han
	Loschbour	LBK	Japanese
	Loschbour	LBK	Papuan
	Loschbour	LBK	Bougainville
	Loschbour	LBK	ASW
	Loschbour	LBK	MXL

We can then run the f3 test like so:

	$ qp3Pop -p f3_par_PaleoNeo > f3_res_PaleoNeo

The output should be familar (except in this case the program helpfully gives us the column headings). 

	qp3Pop: parameter file: f3_par_PaleoNeo
	### THE INPUT PARAMETERS
	##PARAMETER NAME: VALUE
	genotypename: Laz.ASW_MXL.AltaiNea.geno
	snpname: Laz.ASW_MXL.AltaiNea.snp
	indivname: Laz.ASW_MXL.AltaiNea.ind
	popfilename: f3_pop_PaleoNeo
	## qp3Pop version: 412
	packed geno read OK
	end of inpack
	nplist: 17
	number of blocks for block jackknife: 550
	                      Source 1             Source 2               Target           f_3       std. err           Z    SNPs
	 result:             Loschbour                  LBK               Yoruba      0.162301       0.002274      71.361  408343
	 result:             Loschbour                  LBK                Mbuti      0.252859       0.003117      81.133  378650


Which populations appear to show significant evidence of admixture? Are any results borderline, and do such instances make sense?

The third population that is believed to have contributed to modern European genetic ancestry is that of the Bronze Age steppe herders called the Yamnaya. Yamnaya were not in the original Lazaridis et al. 2014 paper, but they were in [Haak et al. Nature. 2015] (https://www.nature.com/articles/nature14317). We have provided you with 9 Yamnaya genomes as PLINK files (***Haak.Yamnaya***). Use the skills you have learnt today to show that the Yamnaya are likely another source population for modern Europeans. If you have time, see if you can examine potential *modern*source populations that may have contributed to African Americans (ASW) and Mexican Americans (MXL).

<!--- 
	$ vi par_convertf_Yamnaya

	genotypename: Haak.Yamnaya.bed
	snpname: Haak.Yamnaya.bim
	indivname: Haak.Yamnaya.fam
	outputformat: PACKEDANCESTRYMAP
	genotypeoutname: Yamnaya.geno
	snpoutname: Yamnaya.snp
	indivoutname: Yamnaya.ind

	$ convertf -p par_convertf_Yamnaya

	$ sed -i 's\???\Yamnaya\g' Yamnaya.ind

	$ vi merge_Yamnaya_par

	geno1: Laz.ASW_MXL.AltaiNea.geno
	snp1:  Laz.ASW_MXL.AltaiNea.snp
	ind1:  Laz.ASW_MXL.AltaiNea.ind
	geno2: Yamnaya.geno
	snp2:  Yamnaya.snp
	ind2:  Yamnaya.ind 
	genooutfilename:   Laz.ASW_MXL.AltaiNea.Yamnaya.geno
	snpoutfilename:    Laz.ASW_MXL.AltaiNea.Yamnaya.snp
	indoutfilename:    Laz.ASW_MXL.AltaiNea.Yamnaya.ind

	$ mergeit -p  merge_Yamnaya_par

	$ vi f3_par_PaleoNeo_Yamnaya

	genotypename: Laz.ASW_MXL.AltaiNea.Yamnaya.geno 
	snpname: Laz.ASW_MXL.AltaiNea.Yamnaya.snp
	indivname: Laz.ASW_MXL.AltaiNea.Yamnaya.ind
	popfilename: f3_pop_PaleoNeo_Yamnaya

	$ vi f3_pop_PaleoNeo_Yamnaya

	Loschbour	Yamnaya	Yoruba
	Loschbour	Yamnaya	Mbuti
	Loschbour	Yamnaya	Sardinian
	Loschbour	Yamnaya	French
	Loschbour	Yamnaya	Orcadian
	Loschbour	Yamnaya	Spanish
	Loschbour	Yamnaya	Sindhi
	Loschbour	Yamnaya	Egyptian
	Loschbour	Yamnaya	Palestinian
	Loschbour	Yamnaya	Karitiana
	Loschbour	Yamnaya	Burusho
	Loschbour	Yamnaya	Han
	Loschbour	Yamnaya	Japanese
	Loschbour	Yamnaya	Papuan
	Loschbour	Yamnaya	Bougainville
	Loschbour	Yamnaya	ASW
	Loschbour	Yamnaya	MXL
	LBK	Yamnaya	Yoruba
	LBK	Yamnaya	Mbuti
	LBK	Yamnaya	Sardinian
	LBK	Yamnaya	French
	LBK	Yamnaya	Orcadian
	LBK	Yamnaya	Spanish
	LBK	Yamnaya	Sindhi
	LBK	Yamnaya	Egyptian
	LBK	Yamnaya	Palestinian
	LBK	Yamnaya	Karitiana
	LBK	Yamnaya	Burusho
	LBK	Yamnaya	Han
	LBK	Yamnaya	Japanese
	LBK	Yamnaya	Papuan
	LBK	Yamnaya	Bougainville
	LBK	Yamnaya	ASW
	LBK	Yamnaya	MXL

	$ qp3Pop -p f3_par_PaleoNeo_Yamnaya > f3_res_PaleoNeo_Yamnaya
--->



##Notes

<sup>1</sup>In this case our input VCF is fairly small (just one sample). It can sometimes be more memory efficient to first convert the VCF file to a transposed PLINK file using '--plink-tped', and then use PLINK to convert from transposed to a regular binary PLINK file.

<sup>2</sup>It is also possible to convert the file to 'EIGENSTRAT' format, but 'PACKEDANCESTRY' map is more compact, saving you valuable hard drive space.

<sup>3</sup>The screen will probably give some warnings starting with "allele funny:". What do you think this is indicating? 

<sup>4</sup>feel free to mess around with this and include your own combinations.

<sup>5</sup>Note that a postive F3 is possible in the case of admixture, but a negative F3 is only possible given admixture, assuming a very simple population history.



