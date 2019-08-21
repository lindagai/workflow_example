################################################################################

#05_phase_whole_chr8_VCF

################################################################################
#0. Setup
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

################################################################################
#NOTE: YOU ONLY NEED TO DO THIS ONCE!

#In Terminal; uncomment out to use

#1. Write shell script to phase the file
cd /users/lgai/latino_datasets/code/R
vim batch.sh

#Type 'i' and copy/paste below line into the file
R CMD BATCH phase.R

#Type `esc`, then type ":wq!" to save

################################################################################

#2. Write shell script to phase the file

#In Terminal; uncomment out to use
vim phase.R

#Type 'i' and copy/paste below line into the file
filepath.beagle4<-"/users/lgai/beagle.r1399.jar"
filepath.vcf<-"/users/lgai/latino_datasets/data/raw_data/latino_chr8_allfiltered.recode.vcf"
filepath.ped<-"/users/lgai/latino_datasets/data/processed_data/peds/latino_peds_cleaned.txt"
filepath.phased.vcf<-"/users/lgai/latino_datasets/data/processed_data/vcfs/PHASED.latino_chr8_allfiltered.recode"

phase.command<-paste0("java -Xmx10000m -jar ", filepath.beagle4,
                      " gt=",filepath.vcf,
                      " ped=",filepath.ped,
                      " out=",filepath.phased.vcf)
phase.command

system(phase.command)

#Type `esc`, then type ":wq!" to save

################################################################################
#3. Submit the job

#In Terminal; uncomment out to use
#Submitted job with this command:
# qsub -cwd -m e -M lindagai@jhu.edu -l mem_free=120G,h_vmem=120G,h_fsize=60G batch.sh

################################################################################