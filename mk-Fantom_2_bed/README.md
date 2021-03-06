# mk-FANTOM_2\_bed - module for converting FANTOM enhancer and promoter data base to VEP friendly format


## About mk-FANTOM_2\_bed

This module converts the FANTOM tab delimited databases into bed files that can be used for [VEP custom annotation](https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html) ([1])

The FANTOM5 Human Enhancers database ([2],[3]) incorporates data produced by the Functional Annotation of the Mammalian Genome consortium (FANTOM). Enhancers were detected using CAP Analysis of Gene Expression. Enhancers show a distinct bidirectional CAGE pattern which allowed its identification in the RNA data produced by FANTOM5 for the majority of cell types and tissues. The database implemented in this pipeline, "enhancer_data_at_2018-05-08_22-04-44.bed", was produced by the slidebase project and it contains coordinates for 32,693 enhancer. It was downloaded from <http://slidebase.binf.ku.dk/human_enhancers/results> on may 22 2018.

The FANTOM5 Human Promoters database ([4],[3]) incorporates data produced by the FANTOM consortium using CAGE sequencing. It contains coordinates for 184,476 promoters from the majority of cell types and tissues. The Database BED file, "promoter_data_at_2018-05-08_22-08-24.bed" was downloaded from <http://slidebase.binf.ku.dk/human_promoters/results> on May 22 2018. 

Coordinates are relative to GRCh37 assembly. The BED file was reformatted to compile with VEP requirements. Annotation files should be stripped of comment lines, sorted in chromosome and position order, compressed using bgzip to be indexed with tabix. The resulting file is then used for annotation.


## Module configuration

The module contains two submodules, one for each kind of regulatory element:

* mk-Enhancer_2\_bed
* mk-Promoter_2\_bed

Within each of them, there are three blocks:

1. Reformat and sort bed file 
2. Compress bed file
3. Index compressed bed with tabix

### Data formats

* Input data
	* Database in BED format
 
* Output data
	* indexed BED file
	
### Software dependencies
* awk, sed, bgzip, tabix.

### Module parameters
* Replacement and substitution of characters (e.g. spaces, field sepatator) with awk and sed
* Selction of relevant fields:
	- chromosome: genomic region
	- start: starting coordinate
	- end: end coordinate
	- enhancerID: ID in the database
	
	or
	
	- promoterID: ID in the database

## mk-Fantom_2\_bed directory structure
````bash
mk-Fantop_2_bed/		##Pipeline main directory
├── mk-Enhancer_2_bed		##Module directory
│   ├── bin		##Executables directory
│   │   └── create-targets		##Script to print files required by this module.
│   ├── data		##Directory with original database file, input file
│   │   └── enhancer_data_at_2018-05-08_22-04-44.bed		##Input file -compressed.
│   ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
│   └── results		##Storage directory for files built by mkfile.  If it does not exist, it is automatically generated by mkfile.
│       ├── enhancer_data_at_2018-05-08_22-04-44.bed		##Intermediate output file.
│       ├── enhancer_data_at_2018-05-08_22-04-44.bed.gz		##output file, sorted and compressed.
│       └── enhancer_data_at_2018-05-08_22-04-44.bed.gz.tbi		##indexed output file.
└── mk-Promoter_2_bed		##Directory for mk sub-module to reformat FANTOM promoter database
    ├── bin		##Executables directory.
    │   └── create-targets		##Script to print every file required by this module.
    ├── data		##Directory with original database file, input file.
    │   └── promoter_data_at_2018-05-08_22-08-24.bed				##Database to be reformatted.
    ├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/create-targets.
    └── results		##Storage directory for files built by mkfile.  If it does not exist, it is automatically generated by mkfile.
        ├── promoter_data_at_2018-05-08_22-08-24.bed		##Intermediate output file.
        ├── promoter_data_at_2018-05-08_22-08-24.bed.gz		##output file, sorted and compressed.
        └── promoter_data_at_2018-05-08_22-08-24.bed.gz.tbi		##indexed output file.


````

## References


[1]: https://www.ncbi.nlm.nih.gov/pubmed/27268795
1. [McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biology. 2016;17:122. doi:10.1186/s13059-016-0974-4.][1]

[2]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5215096/
2. [Andersson R, Gebhard C, Miguel-Escalada I, et al. An atlas of active enhancers across human cell types and tissues. Nature. 2014;507(7493):455-461. doi:10.1038/nature12787.][2]

[3]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5199134/
3. [Ienasescu, H., et al. (2016). On-the-fly selection of cell-specific enhancers, genes, miRNAs and proteins across the human body using SlideBase, Database 2016, doi: 10.1093/database/baw144.][3]

[4]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4529748/
4. [The FANTOM Consortium and the RIKEN PMI and CLST (DGT). A promoter-level mammalian expression atlas. Nature. 2014;507(7493):462-470. doi:10.1038/nature13182.][4]

### Author Info
Developed by [Eugenia Zarza](https://www.researchgate.net/profile/Eugenia_Zarza)(eugenia.zarza@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/). 2018.

