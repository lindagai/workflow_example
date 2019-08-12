########################################
library(VariantAnnotation)
setwd("/users/lgai/latino_datasets/")

########################################

#A. Save common variants VCF
filepath.phased.vcf<-"./data/processed_data/vcfs/PHASED.latino_chr8_allfiltered.recode.vcf"

hg.assembly<-"hg38"

vcf <- VariantAnnotation::readVcf(filepath.phased.vcf, hg.assembly)

#Replace '|' with '/'
geno(vcf)$GT<-gsub("\\|", "\\/",geno(vcf)$GT)

#Replace '1/0' with '0/1'
geno(vcf)$GT[geno(vcf)$GT == "1/0"]<-"0/1"

table(geno(vcf)$GT)

filepath.common.vcf<-"/users/lgai/8q24_project/data/processed_data/vcfs/latino_chr8_allfiltered.for.common.var.vcf"

mv /users/lgai/8q24_project/data/processed_data/vcfs/latino_chr8_allfiltered.for.common.var.vcf /users/lgai/latino_datasets/data/processed_data/vcfs/

writeVcf(vcf,filepath.common.vcf)

########################################

#C. Write bash script

#In Terminal:
# vim batch.sh

#Type 'i', copy/paste the line into below into vim
# R CMD BATCH mend_err_test.R
~
#Type 'esc', and type ':wq' to save/exit

#Run the bash script to submit the job with the following in Terminal

# qsub -cwd -m e -M lindagai@jhu.edu -l mem_free=12G,h_vmem=12G,h_fsize=12G batch.sh

qsub -cwd -m e -M lindagai@jhu.edu -l mem_free=120G,h_vmem=120G,h_fsize=60G batch.sh
