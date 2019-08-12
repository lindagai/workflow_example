################################################################################

#7_examine_mend_test_VCF

################################################################################

ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

################################################################################

#0. Set up
#To run this in terminal from R studio, use option + command + enter

################################################################################

#Packages
library(dplyr)

################################################################################

#1. Examine ggplot

#On laptop console

fp.mend.err.plot <- "/users/lgai/latino_datasets/data/processed_data/chr8_output/chr8.mend.err.plot.png"
filepath.cluster<-"/users/lgai/latino_datasets/data/processed_data/chr8_output/chr8.mend.err.plot.png"
filepath.destination<-" '/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/8q24_example/output' "

scp.command<-paste0("scp lgai@jhpce01.jhsph.edu:", filepath.cluster, " ", filepath.destination)
scp.command
system(scp.command)

################################################################################

#2. Examine Mendelian errors

#On laptop Terminal

#Sign into cluster
ssh -X lgai@jhpce01.jhsph.edu

#Allocate memory
qrsh -l mem_free=15G,h_vmem=16G,h_fsize=18G
module load R
R

fp.mend.err.plot <-"/users/lgai/latino_datasets/data/processed_data/chr8_output/chr8_latino_mend_err.RDS"
mend.errs<-readRDS(fp.mend.err.plot)

mend.errs %>% head

################################################################################

#3. Graph tail end of Mendelian errors if needed

################################################################################

#4. Remove families with Mendelian errors from PED

################################################################################

#5. Remove families with Mendelian errors from VCF