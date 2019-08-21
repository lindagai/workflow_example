################################################################################

#1_rare_var TVTB scratch

################################################################################

#In Terminal; uncomment out to use

# #1. Write shell script to phase the file
# #cd /users/lgai/latino_datasets/code/R/chr8_vcf/common_var_vcf
# cd /users/lgai/latino_datasets/code/R/8q24_vcf/common_var
# vim annot_filter.sh
#
# #Type 'i' and copy/paste below line into the file
# R CMD BATCH trio_geno.R
#
# #Type `esc`, then type ":wq" to save
#
# vim trio_geno.R

################################################################################

#0. Set up
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
#qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
qrsh -l mem_free=30G,h_vmem=30G,h_fsize=30G
# module load R
R

################################################################################

#Packages
library(VariantAnnotation)
library(TVTB)

setwd("/users/lgai/latino_datasets/")

################################################################################

#1. Read in rare variants only - 16 SNPs

fp.ped<-"./data/processed_data/peds/latino_peds_cleaned.txt"
ped <- read.table(fp.ped, header=TRUE, stringsAsFactors = FALSE)
head(ped)

filepath.vcf<-"./data/raw_data/8q24.16snp.08_06_19.vcf"
# filepath.vcf<-"./data/raw_data/8q24.16snp.08_18_19.vcf.bgz"
# filepath.vcf<-"./data/raw_data/8q24.08_18_19.vcf.bgz"
hg.assembly<-"hg38"

vcf <- readVcf(filepath.vcf, hg.assembly)
evcf <- expand(vcf)
evcf

################################################################################

#2. Set up TVTB param

#table(geno(vcf)$GT)
#Can TVTB handle missingness?
# ./.   0/0   0/1   1/1
# 7 10439  1537   737

#Limit the IRange to the 8q24 region
#8q24 is chr8:139864393-141031612

tparam <- TVTBparam(Genotypes(
        ref = "0/0",
        het = c("1/0", "0/1"),
        alt = c("1/1")),
        ranges = GenomicRanges::GRangesList(
                gene.8q24 = GenomicRanges::GRanges(
                        seqnames = "chr8",
                        IRanges::IRanges(
                                start = 139864393, end = 141031612)
                )
        )
)


################################################################################

#3. Define VCF import parameters

#Step 3: Define VCF import parameters
VariantAnnotation::vcfInfo(svp(tparam)) <- vep(tparam)
VariantAnnotation::vcfGeno(svp(tparam)) <- "GT"
svp(tparam)

################################################################################

#4. Define the VCF file to parse

#Step 2: Define the VCF file to parse
vcfFile <- "./data/raw_data/8q24.16snp.08_18_19.vcf.bgz"
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)
tabixVcf

################################################################################

#Skipped these parts

# #Step 1: Import phenotypes
# phenoFile <- system.file(
#         "extdata", "integrated_samples.txt", package = "TVTB")
#
# phenotypes <- S4Vectors::DataFrame(
#         read.table(file = phenoFile, header = TRUE, row.names = 1))
#
# phenotypes

################################################################################

#4.  Import and pre-process variants

#First, must annotate VCF with VEP

#Step 4: Import and pre-process variants

# Import variants as a CollapsedVCF object
vcf <- VariantAnnotation::readVcf(
        tabixVcf, param = tparam, colData = phenotypes)
# Expand into a ExpandedVCF object (bi-allelic records)
vcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)

################################################################################

#5. Add info to the variants

initialInfo <- colnames(info(vcf))
vcf <- addOverallFrequencies(vcf = vcf)
setdiff(colnames(info(vcf)), initialInfo)

################################################################################

#6. Use filter rules to filter

infoR <- VcfInfoRules(
        exprs = list(
                rare = expression(MAF < 0.01 & MAF > 0),
                common = expression(MAF > 0.05)),
        active = c(TRUE, TRUE)
)
infoR

summary(eval(expr = infoR, envir = vcf))

summary(evalSeparately(expr = infoR, envir = vcf))

vepR <- VcfVepRules(exprs = list(
        missense = expression(Consequence %in% c("missense_variant")),
        CADD_gt15 = expression(CADD_PHRED > 15)
))
vepR

################################################################################

#7. Try reading in a file using these rules?
