################################################################################

#3_clean_PED

################################################################################

#0. Set up

ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#Packages
library(dplyr)

setwd("/users/lgai/latino_datasets")

################################################################################

#1. Read in
fp.ped<-"/dcl01/beaty/data/gmkf/latino/june2019/vcfs_filtered/latino_peds.txt"
ped <- read.table(fp.ped,header=FALSE,stringsAsFactors = FALSE)
head(ped)

################################################################################

#2. Correct formatting

colnames(ped) <- c("famid","pid","fatid","motid","sex","affected")

# Sex code ('1' = male, '2' = female)
table(ped$affected)

#Fix coding
ped <- ped %>%
        mutate(sex = ifelse(sex == 1, 1, 0)) %>%
        mutate(affected = ifelse(affected == 3, 1, 0))

head(ped)

################################################################################

#3. Check if VCF IDs match

#Check for IDs:
vcf.pid<-colnames(geno(vcf)$GT)
head(vcf.pid)

head(ped$pid)

sort(setdiff(vcf.pid,ped$pid))
sort(setdiff(ped$pid,vcf.pid))

pids.to.edit <-setdiff(ped$pid,vcf.pid)
pids.to.edit

ped <- ped %>%
        mutate(pid =  ifelse(pid %in% pids.to.edit,
                             paste0(pid,"_2"),pid)) %>%
        mutate(fatid =  ifelse(fatid %in% pids.to.edit,
                               paste0(fatid,"_2"),fatid)) %>%
        mutate(motid =  ifelse(motid %in% pids.to.edit,
                               paste0(motid,"_2"),motid))

################################################################################

#4. Check to ensure no discrepancies between VCF and PED IDs

setdiff(vcf.pid,ped$pid)
setdiff(ped$pid,vcf.pid)
setdiff(ped$fatid,ped$pid)
setdiff(ped$motid,ped$pid)

################################################################################

#5. Write out

filepath.ped.cleaned<-"./data/processed_data/peds/latino_peds_cleaned.txt"
write.table(ped, filepath.ped.cleaned, sep=" ", col.names = TRUE, row.names = FALSE,quote = FALSE)

################################################################################