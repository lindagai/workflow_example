################################################################################

#0_prepare_test_vcf

################################################################################
#0. Setup
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

################################################################################

#1. Prepare your VCF file for read-in

#NOTE: YOU ONLY NEED TO DO THIS ONCE PER VCF!
#In Terminal; uncomment out to use

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G

cd /users/lgai/latino_datasets/data/raw_data

#bgzip the file so you can tabix index it
#this allows you to to use filtering capabilities in VariantAnnotation

module load bcftools
module load tabix
bcftools view 8q24.08_08_19.vcf -Oz -o 8q24.08_18_19.vcf.bgz

tabix -p vcf 8q24.08_18_19.vcf.bgz

################################################################################

# Scratch - 16 snps (no rare var)

################################################################################

#1. Prepare your VCF file for read-in

#NOTE: YOU ONLY NEED TO DO THIS ONCE PER VCF!
#In Terminal; uncomment out to use

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G

cd /users/lgai/latino_datasets/data/raw_data

#bgzip the file so you can tabix index it
#this allows you to to use filtering capabilities in VariantAnnotation

module load bcftools
module load tabix
bcftools view 8q24.16snp.08_06_19.vcf -Oz -o 8q24.16snp.08_18_19.vcf.bgz

tabix -p vcf 8q24.16snp.08_18_19.vcf.bgz

################################################################################
