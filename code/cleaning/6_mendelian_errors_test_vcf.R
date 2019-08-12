################################################################################

#6_mendelian_errors_test_VCF

################################################################################

#Description: Calculates and graphs Mendelian errors for a small test VCF.
# There is no need to run the code here -- it was just used to test an R script
# prior to using it on the complete chr8 VCF.

################################################################################

#0. Set up
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
# ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
# qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
# module load R
# R

################################################################################

#Packages
library(dplyr)
library(trio)
library(VariantAnnotation)

#Set working directory
setwd("/users/lgai/latino_datasets")

################################################################################

#1. Read in files

fp.ped<-"./data/processed_data/peds/latino_peds_cleaned.txt"
ped <- read.table(fp.ped,header=TRUE,stringsAsFactors = FALSE)
head(ped)

filepath.phased.vcf<-"./data/processed_data/vcfs/PHASED.8q24.16snp.08_06_19.vcf.vcf"
hg<-"hg38"
vcf <- readVcf(filepath.phased.vcf, hg)

################################################################################

#2. #Reformat entries into geno format

#Check outcome
table(geno(vcf)$GT)

#Replace '|' with '/'
geno(vcf)$GT<-gsub("\\|", "\\/",geno(vcf)$GT)

#Replace '0/1' with '0/1'
geno(vcf)$GT[geno(vcf)$GT == "1/0"]<-"0/1"

################################################################################

#3. Calculate Mendelian errors

trio.geno<-vcf2geno(vcf,ped)
trio.geno[1:5,1:5]

trio.tmp <- trio::trio.check(dat=trio.geno,is.linkage=FALSE)
trio.tmp$errors %>% head

################################################################################

#4. Graph Mendelian error count and save

if (is.null(trio.tmp$errors)){
        mend.err.out <-trio.tmp$errors
} else {
        #Get Mendelian errors into a dataframe
        mend.err.sorted<-as.data.frame(sort(table(trio.tmp$errors$famid),
                                            decreasing = TRUE))
        colnames(mend.err.sorted)<-c("famid","mend.errors")
        mend.err.out<-mend.err.sorted

        #Graph the Mendelian errors
        ggplot(mend.err.sorted, aes(x=famid,y=mend.errors)) +
                geom_bar(stat="identity") +
                xlab("Family ID") + ylab("Mendelian error count") +
                theme(axis.text.x=element_blank(),
                      axis.ticks.x=element_blank()
                )
        fp.mend.err.plot <- "./data/processed_data/chr8_output/chr8.mend.err.plot.png"
        ggsave()

}

filepath.mend.err<-"./data/processed_data/chr8_output/chr8_latino_mend_err.RDS"
saveRDS(mend.err.out, filepath.mend.err)

################################################################################