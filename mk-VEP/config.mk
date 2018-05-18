##SOFTWARE PATHs	##Should be replaced by direct command calls in the mkfile

VEP="/home/user/ensembl-vep/vep"

##Variant Effect Predictor configurations
#PATH to VEP pluging and databases
#VEPconfig="/home/fperez/.vep/"
VEPconfig="/home/user/.vep" ##Currently the cache is not working. VEP error is: MSG: ERROR: Cannot generate HGVS coordinates (--hgvs) in offline mode without a FASTA file (see --fasta)

#Genome version of the vep cache. GRCh38 or GRCh37 (as long as the vep cache is present)
GENOME_VERSION="GRCh37"

#Database
ANNOTATION_DB="/home/user/Projects/non_coding/non_coding_DB/ORegAnno_2015.hg19.bed.gz"

#Header to describe database annotation
DESCRIPTION="ORegAnno"
#Number of threads to perform annotations
NT="4"
