################################################################################
#0. Set up
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu
#ssh -X lgai@jhpce01.jhsph.edu -o ForwardX11Timeout=336h

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#A. Phase small VCF

filepath.beagle4<-"/users/lgai/beagle.r1399.jar"
filepath.vcf<-"/users/lgai/latino_datasets/data/processed_data/vcfs/8q24.16snp.rarevar.08_19_19.recode.vcf"
filepath.ped<-"/users/lgai/latino_datasets/data/processed_data/peds/latino_peds_cleaned.txt"
filepath.phased.vcf<-"/users/lgai/latino_datasets/data/processed_data/vcfs/PHASED.8q24.16snp.rarevar.08_19_19"

phase.command<-paste0("java -Xmx10000m -jar ", filepath.beagle4,
                      " gt=",filepath.vcf,
                      " ped=",filepath.ped,
                      " out=",filepath.phased.vcf)
phase.command

system(phase.command)

################################################################################
