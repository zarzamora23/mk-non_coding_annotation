# mk-VEP_custom - pipeline for annotation of variants in non-coding regions

## About mk-VEP_costum

Most of the human genome consists of non-protein-coding DNA ([1]). Genome-wide association studies (GWAS) have found links between diseases and single nucleotide polymorphisms (SNPs) occurring in non-coding regions ([2]). Disease-associated variants may affect functional DNA elements that regulate gene expression. It is still challenging to go beyond the statistical associations reported by GWAS and to find the causal relationship between non-coding variants and the affected genes.
Non-coding variants may alter gene expression through multiple mechanims:

* **Promoters**: impose direct impact on transcription initiation and elongation
*  **Enhancers**: Modulation of transcription factor (TF) binding motifs 
*  **Instronic** and **UTR** variants: affect mRNAs altering stability or splicing patterns.
*  **non-coding RNAs**: alter their function or expression

Several consortia have been established to study the effect of non-coding variants through _in silico_ and experimental approaches (e.g. [ENCODE](https://www.encodeproject.org/)). The outcome of these studies has been compiled in large databases. Here we offer a tool to annotate non-coding variants implementing the [Ensemble Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) ([3]) together with four databases:

* [FANTOM5 Human Enhancers database](http://slidebase.binf.ku.dk/human_enhancers/) ([4],[5])
* [FANTOM5 Human Promoters database](http://slidebase.binf.ku.dk/human_promoters/) ([6],[5])
* [Disease Enhancer database](http://biocc.hrbmu.edu.cn/DiseaseEnhancer/) ([7])
* [The Open Regulatory Annotation database (OregAnno)](http://www.oreganno.org/) ([8])

https://doi:10.1093/nar/gkv1203)

## Pipeline configuration

The pipeline is organized in four modules:

### Data formats

* Input data
	* Database in tab delimeted text
	* Variant file in vcf format obtained from whole genome sequencing		
* Output data
	* Database is formatted to bed format with index
	* Annotation produces a tsv and/or vcf file

###Software dependencies
* VEP with libraries to handle bed files

###Module parameters
* See documentation for each module

##mk-VEP-custome directory structure
````
├── mk-DiseaseEnhancer_2_bed
│   ├── bin
│   │   └── create-targets
│   ├── data
│   │   └── enh2disease-1_0_2.txt
│   ├── mkfile
│   └── results
│       ├── enh2disease-1_0_2.bed.gz
│       └── enh2disease-1_0_2.bed.gz.tbi
├── mk-Fantom_2_bed
│   ├── mk-Enhancer_2_bed
│   │   ├── bin
│   │   │   └── create-targets
│   │   ├── data
│   │   │   └── enhancer_data_at_2018-05-08_22-04-44.bed
│   │   ├── mkfile
│   │   └── results
│   │       ├── enhancer_data_at_2018-05-08_22-04-44.bed.gz
│   │       └── enhancer_data_at_2018-05-08_22-04-44.bed.gz.tbi
│   └── mk-Promoter_2_bed
│       ├── bin
│       │   └── create-targets
│       ├── data
│       │   └── promoter_data_at_2018-05-08_22-08-24.bed
│       ├── mkfile
│       └── results
│           ├── promoter_data_at_2018-05-08_22-08-24.bed.gz
│           └── promoter_data_at_2018-05-08_22-08-24.bed.gz.tbi
├── mk-ORegAnno_2_bed
│   ├── bin
│   │   └── create-targets
│   ├── config.mk
│   ├── data
│   │      └── ORegAnno_2015.hg19.tsv
│   ├── mkfile
│   └── results
│       ├── ORegAnno_2015.hg19.bed.gz
│       └── ORegAnno_2015.hg19.bed.gz.tbi
├── mk-VEP
│   ├── config.mk
│   ├── data
│   │   ├── test.vcf
│   │   └── test.vcf.gz
│   ├── mkfile
│   ├── QC
│   │   ├── test_stats.html
│   │   └── test_stats.txt
│   └── results
│       ├── test.vep.ORegAnno.tsv
│       └── test.vep.ORegAnno.vcf
└── README.md
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

