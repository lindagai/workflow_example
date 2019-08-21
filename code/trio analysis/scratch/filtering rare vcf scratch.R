#2. Read in VCF

filepath.vcf<-"./data/raw_data/8q24.16snp.08_06_19.vcf"
# filepath.vcf<-"./data/raw_data/8q24.16snp.08_18_19.vcf.bgz"
# filepath.vcf<-"./data/raw_data/8q24.08_18_19.vcf.bgz"
hg.assembly<-"hg38"




################################################################################

fp.ped<-"./data/processed_data/peds/latino_peds_cleaned.txt"
ped <- read.table(fp.ped, header=TRUE, stringsAsFactors = FALSE)
head(ped)


################################################################################

# Scratch

################################################################################

command.setwd <- "cd /users/lgai/latino_datasets/"
system(command.setwd)
system("pwd")

command.load.vcftools <- "module load vcftools"
system(command.load.vcftools)

filepath.vcf <- "./data/raw_data/8q24.16snp.08_06_19.vcf"
filepath.rare.vcf <- "./data/processed_data/vcfs/8q24.16snp.rarevar.08_19_19.vcf"

command.filter.vcf <- paste0("vcftools --vcf ",
                             filepath.vcf,
                             " --max-maf 0.01",
                             " --non-ref-ac-any 1",
                             " --recode",
                             " --out" ,
                             filepath.rare.vcf)
command.filter.vcf
system(command.filter.vcf)

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
table(geno(vcf)$GT)
# ./.   0/0   0/1   1/1
# 7 10439  1537   737

evcf <- expand(vcf)
evcf

tparam <- TVTBparam(Genotypes(
        ref = "0|0",
        het = c("0|1", "1|0", "0|2", "2|0", "1|2", "2|1"),
        alt = c("1|1", "2|2")),
        ranges = GenomicRanges::GRangesList(
                SLC24A5 = GenomicRanges::GRanges(
                        seqnames = "15",
                        IRanges::IRanges(
                                start = 48413170, end = 48434757)
                )
        )
)

###########

initialInfo <- colnames(info(evcf))
initialInfo
evcf <- addOverallFrequencies(vcf = evcf)
setdiff(colnames(info(evcf)), initialInfo)

###########

#Create VCF filter rule for rare variants
rare.fr <- TVTB::VcfInfoRules(
        exprs = list(
                rare = expression("MAF" < 0.01)
        ),
        active = c(TRUE)
)


rare.fr

rare.filtered <- VariantAnnotation::filterVcf(
        filepath.vcf, hg.assembly, tempfile(), filters=rare.fr
)

rare.filtered

rare.filtered.vcf <- readVcf(rare.filtered, hg.assembly)

rare.filtered.vcf
#dim: 0 795

#??? why isn't it working


################################################################################

#2. Filter rare variant VCF to positions using annotation information

filepath.sm.annovar<-"where/you/want/annovar/to/be.txt"
annovar<-read.table(filepath.sm.annovar,sep="\t",header=TRUE, quote ="")

#Subset to positions of interest
vcf<-vcf[which(start(rowRanges(vcf)) %in% annovar$StartPosition),]

#Write out
filepath.filtered.vcf<-"where/you/want/filtered/VCF/to/be.txt"
writeVCF(rare.filtered, filepath.filtered.vcf)

################################################################################
