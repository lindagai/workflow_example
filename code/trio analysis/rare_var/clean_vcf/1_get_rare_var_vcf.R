################################################################################

#01_prepare_vcf

################################################################################
#0. Setup
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G

################################################################################

#1. Filter VCF to rare variants only (MAF < 0.01) and remove monomorphs
#All code here is to be run in Terminal; uncomment out to use

cd /users/lgai/latino_datasets/
module load vcftools
vcftools --vcf ./data/raw_data/8q24.16snp.08_06_19.vcf --max-maf 0.01 --non-ref-ac-any 1 --recode --out ./data/processed_data/vcfs/8q24.16snp.rarevar.08_19_19

################################################################################
