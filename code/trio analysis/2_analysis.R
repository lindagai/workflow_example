################################################################################

#2_analysis

################################################################################

#In Terminal; uncomment out to use

#1. Write shell script to phase the file
cd /users/lgai/latino_datasets/code/R/chr8_vcf/common_var_vcf
vim trio_analysis.sh

#Type 'i' and copy/paste below line into the file
R CMD BATCH trio_analysis.R

#Type `esc`, then type ":wq" to save

vim trio_analysis.R

# Paste below R code in

################################################################################

#Packages
library(dplyr)
library(VariantAnnotation)
library(trio)

setwd("/users/lgai/latino_datasets/")

################################################################################

#0. Read in

#filepath.trio.geno<-"./data/processed_data/chr8_output/trio.geno.RDS"
#filepath.trio.geno<-"./data/processed_data/8q24_output/8q24.trio.geno.RDS"
filepath.trio.geno<-"./data/processed_data/8q24_output/8q24.16snps.trio.geno.RDS"
trio.geno<-readRDS(filepath.trio.geno)
#dim(trio.geno)
################################################################################

#1. Run gTDT

gTDT.results<-trio::tdt(trio.geno, model = c("additive"))
gTDT.results
dim(trio.geno)

#Create dataframe of genotypic TDT p-values and snp_names
gTDT.df <- data.frame(
        colnames(trio.geno),
        gTDT.results$stat,
        gTDT.results$pval
        )

gTDT.results$pval
gTDT.results$stat

colnames(gTDT.df)<-c("snp","stat","pval")

gTDT.df <- gTDT.df %>%
  arrange(pval) %>%
  mutate(neglogp = -log(pval))

gTDT.df %>% head

#NOTE: I don't know why gTDT only has one pval and one stat...there are 2 SNPs

################################################################################

#2. Run aTDT

aTDT.results<-allelicTDT(trio.geno)

#Create dataframe of allelic TDT p-values and snp_names
aTDT.df<-as.data.frame(cbind(aTDT.results$stat,aTDT.results$pval))
colnames(aTDT.df)<-c("stat","pval")

aTDT.df <- aTDT.df %>%
  arrange(pval) %>%
  mutate(neglogp = -log10(pval))

aTDT.df %>% head

################################################################################

#3. Get position information

#filepath.common.vcf<-"./data/processed_data/vcfs/latino_chr8_allfiltered.for.common.var.vcf"
filepath.common.vcf<-"./data/raw_data/8q24.16snp.08_06_19.vcf"

## Return CHROM, POS, ID and REF from 'fixed' field
# do not read in other fields
hg<-"hg38"
svp <- ScanVcfParam(fixed="NA", info="NA", geno="NA")
vcf.svp <- readVcf(filepath.common.vcf, hg, svp)

snp<-names(vcf.svp)
pos<-start(rowRanges(vcf.svp))
snp.pos.df<-data.frame(snp,pos)
head(snp.pos.df)

gTDT.df<- left_join(gTDT.df,snp.pos.df)
head(gTDT.df)

aTDT.df <- data.frame(gTDT.df[,c("snp","pos")], aTDT.df)
head(aTDT.df)

################################################################################

#4. Save the output

#filepath.gTDT<-"./data/processed_data/chr8_output/gTDT.RDS"
filepath.gTDT<-"./data/processed_data/chr8_output/8q24.16snps.gTDT.RDS"
saveRDS(gTDT.df, filepath.gTDT)

#filepath.aTDT<-"./data/processed_data/chr8_output/aTDT.RDS"
filepath.aTDT<-"./data/processed_data/chr8_output/8q24.16snps.aTDT.RDS"
saveRDS(aTDT.df, filepath.aTDT)

################################################################################

#3. Submit the job

#In Terminal; uncomment out to use
#Submitted job with this command:
# qsub -cwd -m e -M lindagai@jhu.edu -l mem_free=12G,h_vmem=12G,h_fsize=12G trio_analysis.sh

################################################################################
