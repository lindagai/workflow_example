################################################################################

#3_graph_common_var

################################################################################

#0. Set up

ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#Packages
library(ggplot2)
setwd("/users/lgai/latino_datasets")

################################################################################

#1. Read in

fp.gTDT<-"./data/processed_data/chr8_output/8q24.16snps.gTDT.RDS"
fp.aTDT<-"./data/processed_data/chr8_output/8q24.16snps.aTDT.RDS"


gTDT.df <- readRDS(fp.gTDT)
aTDT.df <- readRDS(fp.aTDT)

################################################################################

#2. Plot gTDT

fp.gTDT<-"/users/lgai/latino_datasets/data/processed_data/8q24_output/8q24.16snps.gTDT.pdf"

pdf(fp.gTDT)

ggplot(gTDT.df)+
  geom_point(aes(x=pos,y=neglogp)) +
  labs(title="Genotypic TDT results",
       x="SNP position", y = "-log10p")

dev.off()

################################################################################

#3. Plot aTDT

fp.aTDT<-"/users/lgai/latino_datasets/data/processed_data/8q24_output/8q24.16snps.aTDT.pdf"

pdf(fp.aTDT)

ggplot(aTDT.df)+
  geom_point(aes(x=pos,y=neglogp)) +
  labs(title="Allelic TDT results",
       x="SNP position", y = "-logp")

dev.off()

################################################################################
