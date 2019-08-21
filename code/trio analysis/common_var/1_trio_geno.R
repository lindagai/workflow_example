################################################################################

#1_common_var

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
library(trio)

setwd("/users/lgai/latino_datasets/")

################################################################################

#1. Read in
fp.ped<-"./data/processed_data/peds/latino_peds_cleaned.txt"
ped <- read.table(fp.ped, header=TRUE, stringsAsFactors = FALSE)
head(ped)

# filepath.common.vcf<-"./data/processed_data/vcfs/latino_chr8_allfiltered.for.common.var.vcf"
# filepath.common.vcf<-"./data/raw_data/8q24.08_08_19.vcf"
filepath.common.vcf<-"./data/raw_data/8q24.16snp.08_06_19.vcf"
hg.assembly<-"hg38"
vcf <- readVcf(filepath.common.vcf, hg.assembly)

################################################################################

#2. Get trio matrix

trio.geno <- trio::vcf2geno(vcf,ped)
dim(trio.geno)
trio.geno <- trio::removeSNPs(trio.geno, maf = 0.01)
dim(trio.geno)

#trio.geno[1:5,1:5]
head(trio.geno)

# filepath.trio.geno<-"./data/processed_data/chr8_output/trio.geno.RDS"
filepath.trio.geno<-"./data/processed_data/8q24_output/8q24.16snps.trio.geno.RDS"
saveRDS(trio.geno, filepath.trio.geno)

################################################################################

#3. Submit the job

#In Terminal; uncomment out to use
#Submitted job with this command:
# qsub -cwd -m e -M lindagai@jhu.edu -l mem_free=16G,h_vmem=16G,h_fsize=16G trio_geno.sh

################################################################################
