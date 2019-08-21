################################################################################

#2_annotation_filtering

################################################################################

#0. Set up

#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#Packages
library(dplyr)
library(VariantAnnotation)

setwd("/users/lgai/latino_datasets/")

################################################################################

#1. Read in ANNOVAR report

filepath.annovar<-"/dcl01/beaty/data/gmkf/euro/anno/ANNOVAR_Uploaded/fullannot.gmkf.snvs.c8.T3.ss_ANNOVAR_REPORT.txt"

annovar<-read.table(filepath.annovar, sep="\t", quote ="",
                    header=TRUE, stringsAsFactors=FALSE)

################################################################################

#2. Subset the ANNOVAR to the positions in the VCF

#Get the positions in the VCF by reading in minimal VCF with svp
svp <- ScanVcfParam(fixed="NA", info="NA", geno="NA")
filepath.rare.vcf<-"./data/processed_data/vcfs/8q24.16snp.rarevar.08_19_19.recode.vcf"
vcf.svp <- readVcf(filepath.rare.vcf, "hg38", svp)
snp<-names(vcf.svp)
pos<-start(rowRanges(vcf.svp))

#Complete VCF is large, so let's delete it
rm(vcf.svp)

sm.annovar <- annovar %>% filter(StartPosition %in% pos)

#Complete ANNOVAR is large, so let's delete it
rm(annovar)

################################################################################

#3. Determine which ANNOVAR columns have useful information

#Select useful columns
sm.annovar  <- sm.annovar  %>%
  dplyr::select("StartPosition",
         "EndPosition",
         "ReferenceAllele",
         "AlternativeAllele",
         "Genotype",
         "Quality",
         "TotalDepth",
         "Score_Ljb_pp2hdiv",
         "Score_Ljb_pp2hvar",
         "SiftScore",
         "CADDgt20",
         "WG_GWAVA_score",
         "WG_EIGEN_score"
  )

#Calculate percent-missingness in each annotation score
sm.annovar  %>%
        dplyr::select(
                "Score_Ljb_pp2hdiv",
                "Score_Ljb_pp2hvar",
                "SiftScore",
                "CADDgt20",
                "WG_GWAVA_score",
                "WG_EIGEN_score") %>%
        summarise_each(funs(100*mean(is.na(.))))

# Score_Ljb_pp2hdiv Score_Ljb_pp2hvar SiftScore CADDgt20 WG_GWAVA_score
#               100               100       100      100              0
# WG_EIGEN_score
#               0

#Most values are 100% missing except for GWAVA and EIGEN

#Select functional annotation scores to examine
sm.annovar  %>%
        dplyr::select(
         "Score_Ljb_pp2hdiv",
         "Score_Ljb_pp2hvar",
         "SiftScore",
         "CADDgt20",
         "WG_GWAVA_score",
         "WG_EIGEN_score"
        ) %>%
        lapply(range, na.rm=TRUE)

# $Score_Ljb_pp2hdiv
# [1] 0.001 1.000
#
# $Score_Ljb_pp2hvar
# [1] 0.002 0.997
#
# $SiftScore
# [1] 0.0 0.5
#
# $CADDgt20
# [1] 3.965381 5.221020
#
# $WG_GWAVA_score
# [1] 0.06 0.68
#
# $WG_EIGEN_score
# [1] -1.5687  2.4109

#Select useful columns to save
sm.annovar  <- sm.annovar  %>%
  select("StartPosition",
         "EndPosition",
         "ReferenceAllele",
         "AlternativeAllele",
         "Genotype",
         "Quality",
         "TotalDepth",
         "Score_Ljb_pp2hdiv",
         "Score_Ljb_pp2hvar",
         "SiftScore",
         "CADDgt20",
         "WG_GWAVA_score",
         "WG_EIGEN_score")

################################################################################

#4. Filter ANNOVAR to useful positions, based on recommended score cutoffs from literature

#Select useful columns to save
annovar.filtered <- sm.annovar  %>%
  filter(CADDgt20 > 10 | WG_GWAVA_score > 0.5 | WG_EIGEN_score > 0)

#How many positions are left?
annovar.filtered %>% nrow
# 3217

################################################################################

#5. (OPTIONAL) Subset the ANNOVAR based on results from common variants analysis

#If desired, subset the ANNOVAR to the regions where there is common variant signal
# This looks to be between chr8:129875000-130000000 based on  visual inspection
# But this interval is not in 8q24 proper, and thus not in the VCF that is been read in,
# so we don't use this part of the code in the test VCF

# a<-129875000
# b<-130000000
#
# # 8q24 is chr8:139864393-141031612
# annovar.filtered.peak <- annovar.filtered %>%
#         filter(StartPosition > a & StartPosition < b)
#
# #How many positions are left?
# annovar.filtered.peak %>% nrow

################################################################################

#6. Write out files
fp.annovar.filtered <-"./data/processed_data/annot/ANNOVAR_hg19_filtered.txt"
write.table(annovar.filtered, fp.annovar.filtered, sep="\t",row.names = FALSE,quote = FALSE)
