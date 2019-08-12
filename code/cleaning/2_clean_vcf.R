################################################################################

#2_clean_vcf

################################################################################

#0. Setup

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=16G
R

#Packages
library(VariantAnnotation)

setwd("/users/lgai/latino_datasets")

################################################################################

#1. Read in small portion of VCF as a test file and write it out
fp.vcf<-"./data/raw_data/latino_chr8_allfiltered.recode.vcf.bgz"

#Reading in only a small section (1000 bp range) from 8q24

#8q24 is chr8:139864393-141031612

rng <- GRanges(seqnames="chr8",
               ranges=IRanges(
                       start=c(139864393),
                       end=c(141031612)
               )
)

tab <- TabixFile(fp.vcf)
hg<-"hg38"

#Takes about 5 min
vcf.rng <- readVcf(tab, hg, param=rng)

table(geno(vcf.rng)$GT)
#    ./.     0/0     0/1     1/1
#  14197 9216285  709521  432362

#v. Write out
filepath.sm.vcf<-"./data/raw_data/8q24.08_08_19.vcf"
writeVcf(vcf.rng,filepath.sm.vcf)

################################################################################

#2. Check if it needs to be cleaned

filepath.8q24.vcf<-"./data/raw_data/8q24.08_06_19.vcf"
hg<-"hg38"
vcf <- readVcf(filepath.8q24.vcf, hg)

table(geno(vcf)$GT)
duplicated(start(rowRanges(vcf)))

#The VCF looks clean, so no additional cleaning was performed

################################################################################

#Scratch

rng <- GRanges(seqnames="chr8",
               ranges=IRanges(
                       start=c(129302605),
                       end=c(129303605)
               )
)

tab <- TabixFile(fp.vcf)
hg<-"hg38"

#Takes about 5 min
vcf.rng <- readVcf(tab, hg, param=rng)

#v. Write out
filepath.sm.vcf<-"./data/raw_data/8q24.16snp.08_06_19.vcf"
writeVcf(vcf.rng,filepath.sm.vcf)