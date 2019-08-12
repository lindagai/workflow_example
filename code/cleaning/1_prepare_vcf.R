################################################################################

#01_prepare_vcf

################################################################################
#0. Setup
#To run this in terminal from R studio, use option + command + enter

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

################################################################################

#1. Prepare your VCF file for read-in

#NOTE: YOU ONLY NEED TO DO THIS ONCE PER VCF!
#In Terminal; uncomment out to use

# #note: need to allocate a lot more memory to unzip gzipped file
# qrsh -l mem_free=150G,h_vmem=160G,h_fsize=180G
# cp "/dcl01/beaty/data/gmkf/latino/june2019/vcfs_filtered/latino_chr8_allfiltered.recode.vcf.gz" "/users/lgai/latino_datasets/data/raw_data"

# cd /users/lgai/latino_datasets/data/raw_data
# #bgzip the file so you can tabix index it
# this allows you to to use filtering capabilities in VariantAnnotation
# gunzip latino_chr8_allfiltered.recode.vcf.gz
# module load bcftools
# module load tabix
# bcftools view latino_chr8_allfiltered.recode.vcf -Oz -o latino_chr8_allfiltered.recode.vcf.bgz

# tabix -p vcf latino_chr8_allfiltered.recode.vcf.bgz

################################################################################