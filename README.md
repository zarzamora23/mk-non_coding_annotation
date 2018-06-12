# mk-non-coding-annotation - pipeline for annotation of variants in non-coding regions in addition to coding regions

## About mk-non-coding-annotation

Most of the human genome consists of non-protein-coding DNA ([1]). Genome-wide association studies (GWAS) have found links between diseases and single nucleotide polymorphisms (SNPs) occurring in non-coding regions ([2]). Disease-associated variants may affect functional DNA elements that regulate gene expression. It is still challenging to go beyond the statistical associations reported by GWAS and to find the causal relationship between non-coding variants and the affected genes.
Non-coding variants may alter gene expression through multiple mechanims:

* **Promoters**: impose direct impact on transcription initiation and elongation
* **Enhancers**: Modulation of transcription factor (TF) binding motifs 
* **Intronic** and **UTR** variants: affect mRNAs altering stability or splicing patterns.
* **non-coding RNAs**: alter function or expression

Several consortia have been established to study the effect of non-coding variants through _in silico_ and experimental approaches (e.g. [ENCODE](https://www.encodeproject.org/)). The outcome of these studies has been compiled in large databases. Here we offer a tool to annotate non-coding variants implementing the [Ensemble Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) ([3]) together with four custom databases:

* [FANTOM5 Human Enhancers database](http://slidebase.binf.ku.dk/human_enhancers/results/) ([4],[5])
	Database incorporating data produced by the Functional Annotation of the Mammalian Genome consortium (FANTOM). Enhancers were detected using CAP Analysis of Gene Expression. Enhancers show a distinct bidirectional CAGE pattern which allowed its identification in the RNA data produced by FANTOM5 for the majority of cell types and tissues. The database implemented in this pipeline was produced by the slidebase project and it contains coordinates for 32,693. It was downloaded from (http://slidebase.binf.ku.dk/human_enhancers/results) on may 22 2018. 
	
 * [FANTOM5 Human Promoters database](http://slidebase.binf.ku.dk/human_promoters) ([6],[5])
	This database incorporates data produced by the FANTOM consortium using CAGE sequencing. It contains coordinates for 184,476 promoters from the majority of cell types and tissues. The corresponding BED file was downloaded from (http://slidebase.binf.ku.dk/human_promoters/results) on May 22 2018. In both FANTOM databases, coordinates are relative to GRCh37 assembly.
	
* [Disease Enhancer database](http://biocc.hrbmu.edu.cn/DiseaseEnhancer/) ([7])
	This database provides a comprehensive map of manually curated disease-associated enhancers (1059 disease-associated enhancers in 167 human diseases, involving 896 unique enhancer-gene interactions. Coordinates are relative to the GRCh37 assembly.
	
* [The Open Regulatory Annotation database (OregAnno)](http://www.oreganno.org/) ([8])
	ORegAnno is a crowd-curated regulatory annotation database. It contains information about regulatory regions, transcription factor binding sites, RNA binding sites, regulatory variants, haplotypes, and other regulatory elements. Coordinates in the implemented database are relative to the GRCh37 assembly. However, the original database contains information for other species and for assembly hg18 and GRCh38. The implemented database is the 2015.19.15.tsv version and was downloaded from http://www.oreganno.org/dump/ on May 22 2018 .

## Pipeline configuration

The pipeline is organized in four modules:

* **mk-VEP-custom-multiple**: Creates custom annotation according to custom databases in addition to ensembl data 
* **mk-ORegAnno_2_bed**: Converts ORegAnno database to BED format and indexes it
* **mk-DiseaseEnhancer_2_bed**: Converts DiseaseEnhancer database to BED format and indexes it
* **mk-Fantom_2_bed**: Comprises two submodules, one for each regulatory elements types
	* mk-Enhancer_2_bed: Converts FANTOM5 enhancer database to BED format and indexes it 
	* mk-Promoter_2_bed: Converts FANTOM5 promoter database to BED format and indexes it

### Data formats

* Input data
	* Database in tab delimeted format
	* Chromosome-sorted variant file in vcf or VEP format. It is essential that this file is sorted.
* Output data
	* Database is formatted to BED format with index
	* Annotation produces a tsv and/or vcf file

### Software dependencies
* VEP with Bio::DB::HTS::Tabix module installed to handle bed files, and modules DBI, build DBD with mysql.
* See mk-VEP-custom-multiple for installing instruccions.

### Module parameters
* See documentation for each module

## mk-non-coding-annotation directory structure
````bash
mk-non-custom-annotation-master/		##Pipeline main directory.
├── mk-DiseaseEnhancer_2_bed		##Directory for mk module to reformat DiseaseEnhancer database.
│   ├── bin		##Executables directory.
│   │   └── create-targets		##Script to print every file required by this module.
│   ├── data		##Directory with original database file, input file.
│   │   └── enh2disease-1_0_2.txt		##Database to be reformatted.
│   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/create-targets.
│   ├── README.md	##File explaining how to run this module and expected output.
│   └── results		##Storage directory for files built by mkfile.  If it does not exist, it is automatically generated by mkfile.
│       ├── enh2disease-1_0_2.bed		##Intermediate output file.
│       ├── enh2disease-1_0_2.bed.gz		##output file, sorted and compressed.
│       └── enh2disease-1_0_2.bed.gz.tbi		##indexed output file.
├── mk-Fantom_2_bed		##Directory for mk module to reformat FANTOM databases.
│   ├── mk-Enhancer_2_bed		##Directory for mk sub-module to reformat FANTOM. Enhancer database
│   │   ├── bin		##Executables directory.
│   │   │   └── create-targets		##Script to print every file required by this module.
│   │   ├── data		##Directory with original database file, input file.
│   │   │   └── enhancer_data_at_2018-05-08_22-04-44.bed				##Database to be reformatted.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/create-targets.
│   │   └── results		##Storage directory for files built by mkfile.  If it does not exist, it is automatically generated by mkfile.
│   │       ├── enhancer_data_at_2018-05-08_22-04-44.bed		##Intermediate output file.
│   │       ├── enhancer_data_at_2018-05-08_22-04-44.bed.gz		##output file, sorted and compressed.
│   │       └── enhancer_data_at_2018-05-08_22-04-44.bed.gz.tbi		##indexed output file.
│   ├── mk-Promoter_2_bed		##Directory for mk sub-module to reformat FANTOM promoter database
│   │   ├── bin		##Executables directory.
│   │   │   └── create-targets		##Script to print every file required by this module.
│   │   ├── data		##Directory with original database file, input file.
│   │   │   └── promoter_data_at_2018-05-08_22-08-24.bed				##Database to be reformatted.
│   │   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/create-targets.
│   │   └── results		##Storage directory for files built by mkfile.  If it does not exist, it is automatically generated by mkfile.
│   │       ├── promoter_data_at_2018-05-08_22-08-24.bed		##Intermediate output file.
│   │       ├── promoter_data_at_2018-05-08_22-08-24.bed.gz		##output file, sorted and compressed.
│   │       └── promoter_data_at_2018-05-08_22-08-24.bed.gz.tbi		##indexed output file.
│   └── README.md	##File explaining how to run this module and expected output.
├── mk-ORegAnno_2_bed		##Directory for mk module to reformat ORegAnno database
│   ├── bin		##Executables directory.
│   │   └── create-targets		##Script to print every file required by this module.
│   ├── config.mk		##Configuration file for this module
│   ├── data		##Directory with original database file, input file.
│   │   └── ORegAnno_2015.tsv				##Database to be reformatted.
│   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/create-targets.
│   ├── README.md	##File explaining how to run this module and expected output.
│   └── results		##Storage directory for files built by mkfile.  If it does not exist, it is automatically generated by mkfile.
│       ├── ORegAnno_2015.hg19.bed.gz		##output file, sorted and compressed.
│       └── ORegAnno_2015.hg19.bed.gz.tbi		##indexed output file.
├── mk-VEP-custom-multiple	##Directory for mk module to perform custom annotation and to format output
│   ├── bin		##Executables directory.
│   │   └── create-targets		##Script to print every file required by this module.
│   ├── config.mk		##Configuration file for this module
│   ├── data		##Directory with original database file, input file.
│   │   └── *.vcf.gz		##Input file -compressed.
│   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   ├── QC		##Directory with summary statistics and a list of variants that were not annotated. If it does not exist, it is automatically generated by mkfile.
│   │   ├── *.custom.stats.html		##File in html format showing summary statistics, produced together with the tsv annotated file.
│   │   ├── *.custom.stats.txt		##File in text format showing summary statistics, produced together with the vcf annotated file.
│   │   ├── *.failed_variants.txt	##File containing a list of varinats that were not annotated. The file results from checking for discrepancies between tsv and vcf annotated files.
│   │   ├── *.tsv_variants.txt		##Intermediate file to produce *failed_variants.txt file. Contains list of variants annotated in the tsv file.
│   │   └── *.vcf_variants.txt		##Intermediate file to produce *failed_variants.txt file. Contains list of variants annotated in the vcf file.
│   ├── README.md		##Document describing the pipeline.
│   └── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated.
│       ├── *.columns.tmp.tsv		##Intermediate file. Annotation in tsv fomat without comment headers, just header with column labels.
│       ├── *.DiseaseEnhancer_column.tmp.txt	##Intermediate file. Column containing DiseaseEnhancer annotation only.
│       ├── *.FANTOM_Enhancer_column.tmp.txt	##Intermediate file. Column containing FANTOM Enhancer only.
│       ├── *.FANTOM_Promoter_column.tmp.txt	##Intermediate file. Column containing FANTOM promoter annotation only.
│       ├── *.only_VEP_columns.tmp.tsv	##Intermediate file. Column containing VEP annotation only.
│       ├── *.ORegAnno_column.tmp.txt	##Intermediate file. Column containing ORegAnno annotation only.
│       ├── *.vep.custom_multiple.tsv		##Resulting VEP custom annotation with four databases in tsv format, with the information obtained from a given custom database contained in only one column. 
│       ├── *.vep.custom_multiple.tsv.build_warnings.txt		##File with warning messages generated during the VEP custom annotation with four databases in tsv format.
│       ├── *.vep.custom_multiple.vcf		##Resulting VEP custom annotation with four databases in vcf format.
│       ├── *.vep.custom_multiple.vcf.build_warnings.txt		##File with warning messages generated during the VEP custom annotation with four databases in vcf format.
│       └── *.vep.custom_reformatted.tsv		##Final annotation file in tsv format. Fields in each database are displayed separetely in columns with the corresponding labels. As this file was produced modifying the *.vep.custom_multiple.tsv file, the failed variants are not included and thus, the number of variants may difer from input variant file.
└── README.md		##This document describing the pipeline.

````

## References

[1]: https://www.ncbi.nlm.nih.gov/pubmed/20628352
1. [Alexander RP, Fang G, Rozowsky J, Snyder M, Gerstein MB. Annotating non-coding regions of the genome. Nat Rev Genet. 2010 Aug;11(8):559-71. doi: 10.1038/nrg2814. PubMed PMID: 20628352.][1]

[2]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5529005/
2. [Zhu Y, Tazearslan C, Suh Y. Challenges and progress in interpretation of non-coding genetic variants associated with human disease. Experimental Biology and Medicine. 2017;242(13):1325-1334. doi:10.1177/1535370217713750.][2]

[3]: https://www.ncbi.nlm.nih.gov/pubmed/27268795
3. [McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biology. 2016;17:122. doi:10.1186/s13059-016-0974-4.][3]

[4]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5215096/
4. [Andersson R, Gebhard C, Miguel-Escalada I, et al. An atlas of active enhancers across human cell types and tissues. Nature. 2014;507(7493):455-461. doi:10.1038/nature12787.][4]

[5]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199134/
5. [Ienasescu, H., et al. (2016). On-the-fly selection of cell-specific enhancers, genes, miRNAs and proteins across the human body using SlideBase, Database 2016, doi: 10.1093/database/baw144][5]

[6]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4529748/
6. [The FANTOM Consortium and the RIKEN PMI and CLST (DGT). A promoter-level mammalian expression atlas. Nature. 2014;507(7493):462-470. doi:10.1038/nature13182.][6]

[7]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5753380/
7. [Zhang G, Shi J, Zhu S, et al. DiseaseEnhancer: a resource of human disease-associated enhancer catalog. Nucleic Acids Research. 2018;46(Database issue):D78-D84. doi:10.1093/nar/gkx920.][7]

[8]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4702855/
8. [Lesurf R, Cotto KC, Wang G, et al. ORegAnno 3.0: a community-driven resource for curated regulatory annotation. Nucleic Acids Research. 2016;44(Database issue):D126-D132. doi:10.1093/nar/gkv1203][8]

### Author Info
Developed by [Eugenia Zarza](https://www.researchgate.net/profile/Eugenia_Zarza)(eugenia.zarza@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/). 2018.

