################################################################################

#8_create_rare_and_common_vcfs

################################################################################

ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#Run the next line in Terminal in /users/lgai/latino_datasets/code/R/test_vc
#vim rm.chr8.mend.err.from.vcf.R

#Paste the below code into vim:

########################################

#A. Rare variants VCF - tabix index

#In Terminal

qrsh -l mem_free=150G,h_vmem=160G,h_fsize=60G

gunzip PHASED.latino_chr8_allfiltered.recode.vcf.gz
module load bcftools
bcftools view PHASED.latino_chr8_allfiltered.recode.vcf -Oz PHASED.latino_chr8_allfiltered.recode.vcf.bgz
module load tabix
tabix -p vcf PHASED.latino_chr8_allfiltered.recode.vcf.bgz

#bcftools doesn't work for this

########################################