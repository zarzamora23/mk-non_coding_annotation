##SOFTWARE PATHs	##Should be replaced by direct command calls in the mkfile

VEP="/home/user/ensembl-vep/vep"

##Variant Effect Predictor configurations
#PATH to VEP pluging and databases
VEPconfig="/home/user/.vep" 

#Genome version of the vep cache. GRCh38 or GRCh37 (as long as the vep cache is present)
GENOME_VERSION="GRCh37"

#Databases
DB01="/home/user/Projects/non_coding/non_coding_DB/ORegAnno_2015.hg19.bed.gz"
DB02="/home/user/Projects/non_coding/non_coding_DB/enh2disease-1_0_2.bed.gz"
DB03="/home/user/Projects/non_coding/non_coding_DB/enhancer_data_at_2018-05-08_22-04-44.bed.gz"
DB04="/home/user/Projects/non_coding/non_coding_DB/promoter_data_at_2018-05-08_22-08-24.bed.gz"

#Header to describe database annotation
DESCRIPTION1="ORegAnno"
DESCRIPTION2="DiseaseEnhancer"
DESCRIPTION3="FANTOM_Enhancer"
DESCRIPTION4="FANTOM_Promoter"

