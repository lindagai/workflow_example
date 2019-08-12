################################################################################

#4_phase_test_VCF

################################################################################

#Description: Phases a small test VCF.
# There is no need to run the code here -- it was just used to test an R script
# prior to using it on the complete chr8 VCF.

################################################################################

#0. Set up
ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

setwd("/users/lgai/latino_datasets")

################################################################################

#1. Phase

filepath.beagle4<-"/users/lgai/beagle.r1399.jar"
filepath.vcf<-"./data/raw_data/8q24.16snp.08_06_19.vcf"
filepath.ped<-"./data/processed_data/peds/latino_peds_cleaned.txt"
filepath.phased.vcf<-"./data/processed_data/vcfs/PHASED.8q24.16snp.08_06_19.vcf"

phase.command<-paste0("java -Xmx10000m -jar ", filepath.beagle4,
                      " gt=",filepath.vcf,
                      " ped=",filepath.ped,
                      " out=",filepath.phased.vcf)
phase.command

system(phase.command)

################################################################################