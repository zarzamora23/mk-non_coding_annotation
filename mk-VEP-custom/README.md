# mk-VEP_custom - module to implement custom annotation of non coding variants

## About mk-VEP_custom

The [Ensemble Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) ([1]) determines the effect of variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions. It uses the information compiled by ensembl which is stored in a cache included with the program. Additionally, it allows for the annotation of variants using a custom database.

It require the variants coordinates and the nucleotide changes to find out the:

* Genes and Transcripts affected by the variants
* Location of the variants (e.g. upstream of a transcript, in coding sequence, in non-coding RNA, in regulatory regions)
* Consequence of your variants on the protein sequence (e.g. stop gained, misseense, stop lost, frameshift)
* Known variants that match the input variants, and associated minor allele frequencies from the 1000 Genomes Project
* Information on variants according to a given custom database


## Pipeline configuration

The pipeline is organized in two blocks to produce a vcf or a tab delimited file with ensembl and custom annotation

### Data formats

* Input data
        * vcf file with variants

* Output data
        * vcf file
        * tab delimed text file
        * summary text and html file

### Software dependencies

* VEP can be installed following these steps:

   	1. Download
   	
   		```
         git clone https://github.com/Ensembl/ensembl-vep.git
	```
	
   	2. Install
   	
```
	cd ensembl-vep
	perl INSTALL.pl
```
 	
* VEP with libraries Bio::DB::HTS::Tabix module installed to handle and bed files. These libraries are installed together with VEP

* DBI module. Get the module with cpan

        
        cpan -i DBI
        

* Module Build

       
        cpan Module::Build
       

* Module DBD::mysql otherwise it will run in offline (--offline) mode. Install on Debian ([n]), Ubuntu and derivatives with

        sudo apt-get install libdbd-mysql-perl
   	

[n]: http://search.cpan.org/dist/DBD-mysql/lib/DBD/mysql/INSTALL.pod

* When installing the GHChr37 and 38, VEP may report an error, but has no consequences for the installation:

```
ERROR: For technical reasons this installer is unable to install GRCh37 caches alongside others; please install them
```
## Configuration file

The config.mk file contains important information regarding executables and reference genome. Modify it if needed

    VEP	path to vep executable
    VEPconfig	Variant Effect Predictor configurations
    GENOME_VERSION GRCh37 or GRCH38
    ANNOTATION_DB	Path to the bed formatted database to use for custom annotation
    DESCRIPTION Header for the column showing header annotation

## Module parameters

These are the parameters used to run VEP with a custom database. Information was taken from [VEP's documentation](https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html)

        $VEP	path to executable defined in config.mk file
        -i 	input file
        -o	output file
        -custom	execute custom annotation according to the parameters:
        		$ANNOTATION_DB	defined in config.mk file database used for annotation, 
        		$DESCRIPTION	defined in config.mk filedescription displayed in output, 
        		bed	kind of database file, fotmat type
        		overlap	any annotation that overlaps the variant by even 1bp will be reported
        		0	(default) will output the identifier field if one is found; if none is found, then the coordinates are used instead.
        --stats_file	path stats file
        --cache	use pre-installed cache for annotation
        --offline 	work off line without attempting internet connection
        --dir $VEPconfig	path to VEP directory containing chache
        --assembly genome version to be used, indicated in the config.mk file ($GENOME_VERSION)
        --species study organism, should have genome reference and annotation
        --af	Add the global allele frequency from 1000 Genomes Phase 3 data for any known co-located variant to the output.(AF)
        --af_1kg	Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output. Must be used with --cache.
        --gene_phenotype	Indicates if the overlapped gene is associated with a phenotype, disease or trait.
        --canonical	Adds a flag indicating if the transcript is the canonical transcript for the gene (CANONICAL)
        --check_existing Checks for the existence of known variants that are co-located with your input. By default the alleles are compared and variants on an allele-specific basis - to compare only coordinates, use --no_check_alleles.

        --exclude_null_alleles	Do not include variants with unknown alleles when checking for co-located variants.
        --pick	Pick once line or block of consequence data per variant, including transcript-specific columns. 
        --numbers	Adds affected exon and intron numbering to to output. (EXON,INTRON)
        --regulatory Look for overlaps with regulatory regions. The script can also call if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature (MOTIF_NAME, MOTIF_POS, HIGH_INF_POS, MOTIF_SCORE_CHANGE)
        --biotype	Adds the biotype of the transcript or regulatory feature (BIOTYPE)
        --tab	create output in tab delimited format
        --vcf 	create output in vcf format
        --buffer_size 	Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously. Set this lower to use less memory at the expense of longer run time, and higher to use more memory with a faster run time. Default = 5000


##mk-VEP-custom directory structure
````bash
├── config.mk
├── data
│   ├── test.vcf
│   └── test.vcf.gz
├── mkfile
├── QC
├── Readme.me
└── results

````

## References

[1]: https://www.ncbi.nlm.nih.gov/pubmed/27268795
1. [McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biology. 2016;17:122. doi:10.1186/s1305$



