################################################################################

#2_RV_TDT

################################################################################
#0. Setup
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
R

################################################################################

#Packages
library(VariantAnnotation)
devtools::install_github("lindagai/rvtrio")
library(rvtrio)
install.packages("splitstackshape")
library(splitstackshape)

setwd("/users/lgai/latino_datasets/")

################################################################################

#1. Read in filtered VCF

fp.ped<-"./data/processed_data/peds/latino_peds_cleaned.txt"
ped <- read.table(fp.ped, header=TRUE, stringsAsFactors = FALSE)
head(ped)
dim(ped)

filepath.vcf<-"./data/processed_data/vcfs/PHASED.8q24.16snp.rarevar.08_19_19.vcf.gz"
hg.assembly<-"hg38"
vcf<-VariantAnnotation::readVcf(filepath.vcf, hg.assembly)
dim(geno(vcf)$GT)

################################################################################

#2. Ensure the vcf is coded correctly

table(geno(vcf)$GT)
# 0|0  0|1  1|0  1|1
# 9491   46    2    1

geno(vcf)$GT[geno(vcf)$GT == "0|0"]<-"0/0"
geno(vcf)$GT[geno(vcf)$GT == "0|1"]<-"0/1"
geno(vcf)$GT[geno(vcf)$GT == "1|1"]<-"1/1"
geno(vcf)$GT[geno(vcf)$GT == "1|0"]<-"1/0"

table(geno(vcf)$GT)
# 0/0  0/1  1/0  1/1
# 9491   46    2    1

################################################################################

#3. Run RV_TDT

filepath.to.RV_TDT<-"/users/lgai/rv-tdt/rvTDT"

RV_TDT.results<-rvtrio::RV_TDT(vcf=vcf, vcf.ped = ped, rv.tdt.dir = filepath.to.RV_TDT, upper_cutoff=0.01)
RV_TDT.results





RV_TDT.results<-rvtrio::RV_TDT(vcf=vcf, vcf.ped = ped,
                               rv.tdt.dir = filepath.to.RV_TDT,
                               window.size=5, upper_cutoff=0.01)

#Debugging


################################################################################

#3. Save RV_TDT results

filepath.results<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio_example/example-vcf/RV_TDT.results.360M.windows.txt"
write.table(RV_TDT.results.360M.windows,filepath.results, sep="\t", row.names=FALSE,quote=FALSE)

################################################################################
