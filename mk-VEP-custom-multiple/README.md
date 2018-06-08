# mk-VEP_custom-multiple - module to implement custom annotation of non coding variants

## About mk-VEP_custom

The [Ensemble Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) ([1]) determines the effect of variants (SNPs, insertions, deletions, CNVs or structural variants) on genes, transcripts, and protein sequence, as well as regulatory regions. It uses the information compiled by ensembl which is stored in a cache included with the program. Additionally, it allows for the annotation of variants using custom databases. This pipeline implements annotation with four databases focusing on non-coding regions.

It requires the variants coordinates and the nucleotide changes to find out the:

* Genes and Transcripts affected by the variants
* Location of the variants (e.g. upstream of a transcript, in coding sequence, in non-coding RNA, in regulatory regions)
* Consequence of  variants on the protein sequence (e.g. stop gained, misseense, stop lost, frameshift)
* Known variants that match the input variants, and associated minor allele frequencies from the 1000 Genomes Project
* Information on variants according to four databases that gather information on non-coding regions (i.e. ORegaAnno, Disease Enhancer, FANTOM Enhancer, FANTOM promoter;  see below for details)

## Pipeline configuration

The pipeline is organized in nine  blocks,  to produce annotation in vcf and tab delimited format (2 blocks) , to reformat the resulting tsv file to obtain a spreadsheet freindly file (6 blocks), and one block to find the variants that failed and that were not annotated in the vcf.

### Data formats

* Input data
        -  Chromosome-sorted vcf or VEP format file with variants

* Output data
        - vcf file
        - tab delimed text file
        - summary text and html file
        - file with list of failed variants, see below for details
        
	**WARNING**: Analysing an input file may produce different outcomes in the vcf and tab files when applting the arguments "--dont_skip" and " --allow _non_variant". The vcf wil contain the same variants as in the input file. Failed variants will be listed but will not have annotation. On the contrary, the tsv excludes variants that failed validation, thus a smaller number of variants will be reported. Example of failed variant: 1       155194980       155194980       T/T. 

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
 	
* VEP with libraries Bio::DB::HTS::Tabix module installed to handle bed files. These libraries are installed together with VEP

* DBI module. Get the module with cpan

        
        cpan -i DBI
        

* Module Build

       
        cpan Module::Build
       

* Module DBD::mysql otherwise it will run only in offline (--offline) mode. [Install](http://search.cpan.org/dist/DBD-mysql/lib/DBD/mysql/INSTALL.pod) on Debian, Ubuntu and derivatives with

        sudo apt-get install libdbd-mysql-perl
   	
* When installing the GHChr37 and 38, VEP may report an error, but has no consequences for the installation:

```
ERROR: For technical reasons this installer is unable to install GRCh37 caches alongside other
```
## Configuration file

The config.mk file contains important information regarding loaction of executables, reference genome version to be used, and database. Modify it if needed

    VEP	path to vep executable
    VEPconfig	Variant Effect Predictor configurations, indicate path if not installed in default location
    GENOME_VERSION GRCh37 or GRCH38. Will search for them in VEPconfig
    DB01-DB04	Paths to BED formatted databases to use for custom annotation.
    	
    	DB01="/path/to/database/ORegAnno_2015.hg19.bed.gz"
    	DB02="/path/to/database/enh2disease-1_0_2.bed.gz"
    	DB03="/path/to/database/enhancer_data_at_2018-05-08_22-04-44.bed.gz"
    	DB04="/path/to/database/promoter_data_at_2018-05-08_22-08-24.bed.gz"

    DESCRIPTION1 - DESCRIPTION4 Headers for column showing annotation
    
		DESCRIPTION1="ORegAnno"
		DESCRIPTION2="DiseaseEnhancer"
		DESCRIPTION3="FANTOM_Enhancer"
		DESCRIPTION4="FANTOM_Promoter"


## Module parameters

These are the parameters used to run VEP with a custom database. Information was taken from [VEP's documentation](https:/www.ensembl.org/info/docs/tools/vep/script/vep_options.html). The resulting column header is indicated in parenthesis.

        $VEP \ -> Path to executable defined in config.mk file
        -i  \ -> input file (in vcf or VEP format, must be sorted by chromosome)
        -o  \ -> output file
        -custom  \ -> execute custom annotation according to the parameters:
        		$DB_0X  \ -> defined in config.mk file, database used for annotation. Currently, it ould range from DB_01 - DB_04 as we have implemented four databases
        		$DESCRIPTION  \ -> defined in config.mk file, annotation description that will be displayed in output, 
        		BED  \ -> kind of database file to be used for annotation, fotmat type; GFF, GTF, VCF and bigWig are also supported,
        		overlap  \ -> any annotation that overlaps the variant by even 1bp will be reported; alternative parameter is "exact" where only annotations whose coordinated match exactly those of the variant will be reported
        		0 \ -> (default) will output the identifier field if one is found; if none is found, then the coordinates are used instead; alternattive parameter is 1, when used VEP will output the coordinates of an overlapping feature.
        --stats_file \ -> path to output stats file
        --cache \ -> use pre-installed cache for annotation
        --offline 	work off line without attempting internet connection
        --dir $VEPconfig \ -> path to VEP directory containing cache
        --assembly \ -> genome assembly version to be used, indicated in the config.mk file as the variable $GENOME_VERSION
        --species \ -> study organism, should have genome reference and annotation
        --af \ -> Add the global allele frequency from [1000 Genomes Phase 3](http://www.internationalgenome.org/ data) for any known co-located variant to the output.(AF)
        --af_1kg \ -> Add allele frequency from continental populations (AFR - Africa-, AMR - America -,EAS -East Asia- ,EUR -Europe-,SAS -South Asia) of 1000 Genomes Phase 3 to the output. Must be used with --cache.
        --gene_phenotype \ -> Indicates if the overlapped gene is associated with a phenotype, disease or trait.
        --canonical \ -> Adds a flag indicating if the transcript is the canonical transcript for the gene (CANONICAL)
        --check_existing \ -> Checks for the existence of known variants that are co-located with your input. By default the alleles are compared and variants on an allele-specific basis - to compare only coordinates, use --no_check_alleles.
        --exclude_null_alleles \ -> Do not include variants with unknown alleles when checking for co-located variants.
        --pick \ -> Pick once line or block of consequence data per variant, including transcript-specific columns. 
        --numbers \ -> Adds affected exon and intron numbering to to output. (EXON,INTRON)
        --regulatory \ -> Look for overlaps with regulatory regions in the Ensembl database. The script can also call if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature (MOTIF_NAME, MOTIF_POS, HIGH_INF_POS, MOTIF_SCORE_CHANGE)
        --biotype \ -> Adds the biotype of the transcript or regulatory feature (BIOTYPE)
        --variant_class \ -> Output the Sequence Ontology variant class (VARIANT_CLASS)
        --tab \ -> create output in tab delimited format
        #or
        --vcf  \ -> create output in vcf format
        --sift [p|s|b] \ -> SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both.  (SIFT)
        --polyphen [p|s|b] \ -> predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both. VEP uses the humVar score by default - use --humdiv to retrieve the humDiv score. 
        --dont_skip \ -> Do not skip input variants that fail validation, e.g. those that fall on unrecognised sequences 	 
	--allow_non_variant \ -> lines where the ALT allele is null, are printed in the VCF output with no consequence data added. 
        --buffer_size \ -> Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously. Set this lower to use less memory at the expense of longer run time, and higher to use more memory with a faster run time. Default = 5000

## Information provided in the output tsv file per column

* **Uploaded_variation**: variant ID in the format "chromosome\_start\_alleles", based on variant coordinates
* **Location**: in standard coordinate format (chr:start or chr:start-end)
* **Allele**: the variant allele used to calculate the consequence (variant in the sample)
* **Gene**: Ensembl stable ID of affected gene
* **Feature**: Ensembl stable ID of feature
* **Feature_type**: type of feature (i.e. Transcript, RegulatoryFeature, MotifFeature)
* **[Consequence](https://www.ensembl.org/info/genome/variation/predicted_data.html#consequences)**: VEP identifies each Ensembl transcript overlapping a mapped variant, and then predicts the effect that each allele may have on the transcript. VEP employs terms defined by the [Sequence Onthology](http://www.sequenceontology.org/) to describe the features and attributes of biological sequence. We show below the meaning of each of the possible outcomes and their impact in the protein (see next point for details). Their location in respect to the transcript is shown in Fig. 1.
	 
	- Transcript ablation: A feature loss whereby the deleted region includes a transcript feature, high impact (HI)
	- Splice acceptor variant: A splice variant that changes the 2 base region at the 3' end of an intron, HI
	- Splice donor variant: A splice variant that changes the 2 base region at the 5' end of an intron, HI
	- Stop gained: A sequence variant whereby at least one base of a codon is changed, resulting in a premature stop codon, leading to a shortened transcript, HI
	- Frameshift variant: A sequence variant which causes a disruption of the translational reading frame, because the number of nucleotides inserted or deleted is not a multiple of three, HI
	- Stop lost: A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript, HI
	- Start lost: A codon variant that changes at least one base of the canonical start codon, HI
	- Transcript amplification: A feature amplification of a region containing a transcript, HI
	- Inframe insertion: An inframe non synonymous variant that inserts bases into in the coding sequence, moderate impact (MI)
	- Inframe deletion: An inframe non synonymous variant that deletes bases from the coding sequence, MI
	- Missense variant: A sequence variant, that changes one or more bases, resulting in a different amino acid sequence but where the length is preserved, MI
	- Protein altering variant: A sequence_variant which is predicted to change the protein encoded in the coding sequence, MI
	- Splice region variant: A sequence variant in which a change has occurred within the region of the splice site, either within 1-3 bases of the exon or 3-8 bases of the intron, low impact (LI)
	- Incomplete terminal codon variant: A sequence variant where at least one base of the final codon of an incompletely annotated transcript is changed, LI
	- Start retained variant: A sequence variant where at least one base in the start codon is changed, but the start remains, LI
	- Stop retained variant: A sequence variant where at least one base in the terminator codon is changed, but the terminator remains, LI
	- Synonymous variant: A sequence variant where there is no resulting change to the encoded amino acid, LI
	- Coding sequence variant: A sequence variant that changes the coding sequence, modifier (M)
	- Mature miRNA variant: A transcript variant located with the sequence of the mature miRNA, M
	- 5 prime UTR variant: A UTR variant of the 5' UTR, M
	- 3 prime UTR variant: A UTR variant of the 3' UTR, M
	- Non coding transcript exon variant: A sequence variant that changes non-coding exon sequence in a non-coding transcript, M
	- Intron variant: A transcript variant occurring within an intron, M
	- NMD transcript variant: A variant in a transcript that is the target of NMD, M
	- Non coding transcript variant: A transcript variant of a non coding RNA gene, M
	- Upstream gene variant: A sequence variant located 5' of a gene, M
	- Downstream gene variant: A sequence variant located 3' of a gene, M
	- TFBS ablation: A feature ablation whereby the deleted region includes a transcription factor binding site, M
	- TFBS amplification: A feature amplification of a region containing a transcription factor binding site,	M
	- TF binding site variant: A sequence variant located within a transcription factor binding site, M
	- Regulatory region ablation: A feature ablation whereby the deleted region includes a regulatory region, M
	- Regulatory region amplification: A feature amplification of a region containing a regulatory region, M
	- Feature elongation: A sequence variant that causes the extension of a genomic feature, with regard to the reference sequence, M
	- Regulatory region variant: A sequence variant located within a regulatory region, M
	- Feature truncation: A sequence variant that causes the reduction of a genomic feature, with regard to the reference sequence, M
	- Intergenic variant: A sequence variant located in the intergenic region, between genes, M
	
* **[IMPACT](https://www.ensembl.org/Help/Glossary?id=535)**: the impact modifier for the consequence type, according to the scale suggested by VEP:
	- High impact: Disruptive impact in the protein, may cause protein truncation, loss of function or triggering nonsense mediated decay
	- Moderate impact: It might change protein effectiveness
	- Low impact: Mostly harmless, unlikely to change protein behaviour
	- Modifier: Usually non-coding variants or variants affecting non-coding genes, where predictions are difficult or there is no evidence of impact.
	
* **cDNA_position**: relative position of base pair in cDNA sequence
* **CDS_position**: relative position of base pair in coding sequence
* **Protein_position**: relative position of amino acid in protein
* **Amino_acids**: only given if the variant affects the protein-coding sequence
* **Codons**: the alternative codons with the variant base in upper case
* **Existing_variation**: known identifier of existing variant
* **DISTANCE**: Shortest distance from variant to transcript
* **STRAND**: the DNA strand (1 or -1) on which the transcript/feature lies
* **FLAGS**: transcript quality flags: cds_start_NF: CDS 5' incomplete; cds_end_NF: CDS 3' incomplete
* **BIOTYPE**: Biotype of transcript or regulatory feature
* **VARIANT_CLASS**: [Sequence Ontology variant class](https://www.ensembl.org/info/genome/variation/data_description.html#classes). We show below the description of the variant class annotation, first according to the Sequence Ontology classification and in parenthesis as the term used by VEP. The default is to use the Sequence Ontology.

	- SNV: SNVs are single nucleotide positions in genomic DNA at which different sequence alternatives exist. (Variant)
	- genetic_marker: A measurable sequence feature that varies within a population. (Variant)
	- substitution: A sequence alteration where the length of the change in the variant is the same as that of the reference. (Variant)
	- tandem_repeat: Two or more adjacent copies of a region (of length greater than 1). (Variant)
	- Alu_insertion: An insertion of sequence from the Alu family of mobile elements. (SV)
	- complex_structural_alteration: A structural sequence alteration or rearrangement encompassing one or more genome fragments, with 4 or more breakpoints. (SV)
	- complex_substitution: When no simple or well defined DNA mutation event describes the observed DNA change, the keyword \"complex\" should be used. Usually there are multiple equally plausible explanations for the change. (SV)
	- copy_number_gain: A sequence alteration whereby the copy number of a given region is greater than the reference sequence. (SV)
	- copy_number_loss: A sequence alteration whereby the copy number of a given region is less than the reference sequence. (SV)
	- copy_number_variation: A variation that increases or decreases the copy number of a given region. (SV)
	- duplication: An insertion which derives from, or is identical in sequence to, nucleotides present at a known location in the genome.	SV
	- interchromosomal_breakpoint: A rearrangement breakpoint between two different chromosomes. (SV)
	- interchromosomal_translocation: A translocation where the regions involved are from different chromosomes. (SV)
	- intrachromosomal_breakpoint: A rearrangement breakpoint within the same chromosome. (SV)
	- intrachromosomal_translocation: A translocation where the regions involved are from the same chromosome. (SV)
	- inversion: A continuous nucleotide sequence is inverted in the same position. (SV)
	- loss_of_heterozygosity: A functional variant whereby the sequence alteration causes a loss of function of one allele of a gene. (SV)
	- mobile_element_deletion: A deletion of a mobile element when comparing a reference sequence (has mobile element) to a individual sequence (does not have mobile element). (SV)
	- mobile_element_insertion: A kind of insertion where the inserted sequence is a mobile element. (SV)
	- novel_sequence_insertion: An insertion the sequence of which cannot be mapped to the reference genome. (SV)
	- short_tandem_repeat_variation: A variation that expands or contracts a tandem repeat with regard to a reference. (SV)
	- tandem_duplication: A duplication consisting of 2 identical adjacent regions. (SV)
	- translocation: A region of nucleotide sequence that has translocated to a new position. The observed adjacency of two previously separated regions. (SV)
	- deletion: The point at which one or more contiguous nucleotides were excised.	(Variant, SV)
	- indel: A sequence alteration which included an insertion and a deletion, affecting 2 or more bases. (Variant, SV)
	- insertion: The sequence of one or more nucleotides added between two adjacent nucleotides in the sequence. (Variant, SV)
	- sequence_alteration: A sequence_alteration is a sequence_feature whose extent is the deviation from another sequence. (Variant, SV)
	- probe: A DNA sequence used experimentally to detect the presence or absence of a complementary nucleic acid. (CNV probe)

* **CANONICAL**: a flag indicating if the transcript is denoted as the canonical transcript for this gene
* **SOURCE**: Annotation source in custom annotation. Only printed when using GFF or GTF files.
* **GENE_PHENO**: Indicates if overlapped gene is associated with a phenotype, disease or trait
* **SIFT**: the SIFT prediction and/or score, with both given as prediction(score)
* **PolyPhen**: the PolyPhen prediction and/or score
* **INTRON**: the intron number (out of total number)
* **EXON**: the exon number (out of total number)
* **AF**: Frequency of existing variant in 1000 Genomes
* **AFR_AF**: Frequency of existing variant in 1000 Genomes combined African population
* **AMR_AF**: Frequency of existing variant in 1000 Genomes combined American population
* **EAS_AF**: Frequency of existing variant in 1000 Genomes combined East Asian population
* **EUR_AF**: Frequency of existing variant in 1000 Genomes combined European population
* **SAS_AF**: Frequency of existing variant in 1000 Genomes combined South Asian population
* **CLIN_SIG**: Clinical significance of the dbSNP variant according to the ClinVar database (association, beningn, confers sensitivity, drug response, likely benign, likely pathogenic, not provided, other, pathogenic, protective, risk factor, uncertain significance) 
* **SOMATIC**: Somatic status of existing variant(s); multiple values correspond to multiple values in the Existing_variation field
* **PHENO**: Indicates if existing variant is associated with a phenotype, disease or trait; multiple values correspond to multiple values in the Existing_variation field
* **MOTIF_NAME**: the source and identifier of a transcription factor binding profile aligned at this position
* **MOTIF_POS**: The relative position of the variation in the aligned TFBP
* **HIGH_INF_POS**: a flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)
* **MOTIF_SCORE_CHANGE**: The difference in motif score of the reference and variant sequences for the TFBP
* **OREGgene_symbol**: Gene symbol in the OregAnno database
* **OREGstrand**: Strand where the variant occurs according to OregAnno database
* **OREGtype**: Regulatory type (e.g Transcription factor binding site)
* **OREGexperimental_outcome**: Each record is associated to a positive, neutral or negative outcome based on the experimental results from the primary reference.
* **OREGdatabase_ID**: stable OregAnno identifier
* **OREGNCBI**: Pubmed accession for publication supporting variant
* **DEgene_symbol**: Gene symbol in the Disease Enhancer database
* **DEdisease**: Associated disease
* **DEdatabase_ID**: Disease Enhancer data base ID
* **FANTOM_Enhancer**: Enhancer coordinates according to FANTOM database
* **FANTOM_Promoter**: Promoter coordinates according to FANTOM database

![Consequences](/home/user/Documentos/ReadMes/consequences.png  "Location of consequence terms")
Fig. 1 . Location of VEP consequence terms relative to transcript. Taken from https://www.ensembl.org/info/genome/variation/predicted_data.html#consequences 

## mk-VEP-custom-multiple
````bash
mk-VEP-custom-multiple/		##Pipeline main directory.
├── bin		##Executables directory
│   └── create-targets		##Script to print files required by this module.
├── config.mk		##Configuration file for this module.
├── data		##Directory containing input file
│   └── *.vcf.gz		##Input file -compressed.
├── mkfile		##File in mk format, specifying the rules for building every result requested by bin/targets.
├── QC		##Directory with summary statistics and a list of variants that were not annotated. If it does not exist, it is automatically generated by mkfile.
│   ├── *.custom.stats.html		##File in html format showing summary statistics, produced together with the tsv annotated file.
│   ├── *.custom.stats.txt		##File in text format showing summary statistics, produced together with the vcf annotated file.
│   ├── *.failed_variants.txt	##File containing a list of varinats that were not annotated. The file results from checking for discrepancies between tsv and vcf annotated files.
│   ├── *.tsv_variants.txt		##Intermediate file to produce *failed_variants.txt file. Contains list of variants annotated in the tsv file.
│   └── *.vcf_variants.txt		##Intermediate file to produce *failed_variants.txt file. Contains list of variants annotated in the vcf file.
├── README.md		##This document describing the pipeline.
└── results		##Storage directory for files built by mkfile. If it does not exist, it is automatically generated.
    ├── *.columns.tmp.tsv		##Intermediate file. Annotation in tsv fomat without comment headers, just header with column labels.
    ├── *.DiseaseEnhancer_column.tmp.txt	##Intermediate file. Column containing DiseaseEnhancer annotation only.
    ├── *.FANTOM_Enhancer_column.tmp.txt	##Intermediate file. Column containing FANTOM Enhancer only.
    ├── *.FANTOM_Promoter_column.tmp.txt	##Intermediate file. Column containing FANTOM promoter annotation only.
    ├── *.only_VEP_columns.tmp.tsv	##Intermediate file. Column containing VEP annotation only.
    ├── *.ORegAnno_column.tmp.txt	##Intermediate file. Column containing ORegAnno annotation only.
    ├── *.vep.custom_multiple.tsv		##Resulting VEP custom annotation with four databases in tsv format, with the information obtained from a given custom database contained in only one column. 
    ├── *.vep.custom_multiple.tsv.build_warnings.txt		##File with warning messages generated during the VEP custom annotation with four databases in tsv format.
    ├── *.vep.custom_multiple.vcf		##Resulting VEP custom annotation with four databases in vcf format
    ├── *.vep.custom_multiple.vcf.build_warnings.txt		##File with warning messages generated during the VEP custom annotation with four databases in vcf format.
    └── *.vep.custom_reformatted.tsv		##Final annotation file in tsv format. Fields in each database are displayed separetely in columns with the corresponding labels. As this file was produced modifying the *.vep.custom_multiple.tsv file, the failed variants are not included and thus, the number of variants may difer from input variant file.
````

## References

[1]: https://www.ncbi.nlm.nih.gov/pubmed/27268795
1.  [McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biology. 2016;17:122. doi:10.1186/s13059-016-0974-4.][1]

### Author Info
Developed by [Eugenia Zarza](https://www.researchgate.net/profile/Eugenia_Zarza)(eugenia.zarza@gmail.com) for [Winter Genomics](http://www.wintergenomics.com/). 2018.
