#Test 3: Check input file dimension and entries

########################################################

setwd("/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets")

########################################################

library(VariantAnnotation)
library(rvtrio)
library(splitstackshape)

################################################################################

#1. Check dimensions of Latino dataset - OK

fp.map <- "/Users/lindagai 1/Documents/classes/4th year/Research/window1.1-12M.map"
map <- read.table(fp.map, header=FALSE, stringsAsFactors = FALSE)
map

fp.tped <- "/Users/lindagai 1/Documents/classes/4th year/Research/window1.1-12M.tped"
tped <- read.table(fp.tped, header=FALSE, stringsAsFactors = FALSE)
tped[1:5,1:5]

fp.ped.rvtdt <- "/Users/lindagai 1/Documents/classes/4th year/Research/pedfile.ped"
ped.rvtdt <- read.table(fp.ped.rvtdt, header=FALSE, stringsAsFactors = FALSE)
ped.rvtdt

#Checks
#1. Are the dimensions correct?
dim(map)
#12  3

dim(tped)
#12 1591
795*2 +1 #1591

dim(ped.rvtdt)
#795   6

dim(geno(vcf)$GT)
#12 795

################################################################################

#2. Do the SNP names in tped match the names in the map file? - OK

identical(tped[,1], map[,2])
#TRUE

################################################################################

#3. Are the table values all phased? - OK

filepath.vcf<-"/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/PHASED.8q24.16snp.rarevar.08_19_19.vcf.gz"
hg.assembly<-"hg38"
vcf<-VariantAnnotation::readVcf(filepath.vcf, hg.assembly)
dim(geno(vcf)$GT)

table(geno(vcf)$GT)

geno(vcf)$GT[geno(vcf)$GT == "0|0"]<-"0/0"
geno(vcf)$GT[geno(vcf)$GT == "0|1"]<-"0/1"
geno(vcf)$GT[geno(vcf)$GT == "1|1"]<-"1/1"
geno(vcf)$GT[geno(vcf)$GT == "1|0"]<-"1/0"

table(geno(vcf)$GT)
# 0/0  0/1  1/0  1/1 
9491   46    2    1 

################################################################################

#4. Do the PED and VCF Geno look correct? - OK

head(ped)
# 
# pid      famid      fatid      motid sex affected
# 1 GMKFC130701201 GMKFC15913          0          0   0        0
# 2 GMKFC130700286 GMKFC15927 GMKFC15929 GMKFC15928   0        1
# 3 GMKFC130700968 GMKFC15934          0          0   1        0
# 4 GMKFC130702661 GMKFC15973          0          0   1        0
# 5 GMKFC130704257 GMKFC16007          0          0   0        0
# 6 GMKFC130701320 GMKFC16170 GMKFC16171 GMKFC16172   0        1

geno(vcf)$GT[1:5,1:5]

# GMKFC15209 GMKFC15212 GMKFC15214 GMKFC15257 GMKFC15258
# chr8:129302656_T/C "0/0"      "0/0"      "0/0"      "0/0"      "0/0"     
# rs144566303        "0/0"      "0/0"      "0/0"      "0/0"      "0/0"     
# rs139161223        "0/0"      "0/0"      "0/0"      "0/0"      "0/0"     
# chr8:129302863_G/A "0/0"      "0/0"      "0/0"      "0/0"      "0/0"     
# rs368850832        "0/0"      "0/0"      "0/0"      "0/0"      "0/0" 
