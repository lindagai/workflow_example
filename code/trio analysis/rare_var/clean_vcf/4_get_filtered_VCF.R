################################################################################

#2_RV_TDT

################################################################################

#Packages
library(VariantAnnotation)

setwd("/users/lgai/latino_datasets/")

################################################################################

#1. Read in
fp.annovar.filtered <-"./data/processed_data/annot/ANNOVAR_hg19_filtered.txt"
annovar<-read.table(fp.annovar.filtered, sep="\t",header=TRUE, quote ="")

filepath.phased.vcf<-"./data/processed_data/vcfs/PHASED.8q24.16snp.rarevar.08_19_19.vcf.gz"
hg.assembly<-"hg38"
vcf <- readVcf(filepath.phased.vcf, hg.assembly)

################################################################################

#2. Filter to positions in ANNOVAR

#None of the ANNOVAR positions were in the 16 SNP file, so this code was not run

vcf.filtered<-vcf[which(start(rowRanges(vcf)) %in% annovar$StartPosition),]
dim(geno(vcf.filtered)$GT)

filepath.vcf.filtered<-"./data/processed_data/vcfs/8q24.08_08_19.vcf"
writeVcf(vcf.filtered,filepath.vcf.filtered)
