################################################################################

#2_RV_TDT

################################################################################

#In Terminal; uncomment out to use

#1. Write shell script to phase the file
#cd /users/lgai/latino_datasets/code/R/chr8_vcf/common_var_vcf
cd /users/lgai/latino_datasets/code/R/8q24_vcf/common_var
vim trio_geno.sh

#Type 'i' and copy/paste below line into the file
R CMD BATCH trio_geno.R

#Type `esc`, then type ":wq" to save

vim trio_geno.R

################################################################################

#Packages
library(VariantAnnotation)
library(TVTB)

setwd("/users/lgai/latino_datasets/")

################################################################################

#3. Read in filtered VCF
fp.ped<-"./data/processed_data/peds/latino_peds_cleaned.txt"
ped <- read.table(fp.ped, header=TRUE, stringsAsFactors = FALSE)
head(ped)

filepath.vcf<-"./data/raw_data/8q24.16snp.08_06_19.vcf"
hg.assembly<-"hg38"

################################################################################
