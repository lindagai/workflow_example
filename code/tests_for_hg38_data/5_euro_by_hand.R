#Test 5: Do you still get the error in the Latino dataset if you make the TPED/PED/MAP by hand?

########################################################

setwd("/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets")

########################################################

library(VariantAnnotation)
library(rvtrio)
library(splitstackshape)

########################################################

#1. Read in PED

filepath.vcf.ped<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.phen"

ped<-read.table(filepath.vcf.ped,header=FALSE)
colnames(ped)<-c("pid","famid","fatid","motid","sex","affected")
head(ped)

########################################################

#2. Read in VCF

filepath.vcf<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/data/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"
hg.assembly<-"hg19"
vcf<-VariantAnnotation::readVcf(filepath.vcf, hg.assembly)

########################################################

#3. Check VCF entries
table(geno(vcf)$GT)
# 0/0    0/1    1/1 
# 357734   2290    984 

########################################################

#4. Make the MAP, TPED, and PED file by hand

snp.pos.df<-as.data.frame(cbind(names(vcf),start(rowRanges(vcf))))
colnames(snp.pos.df)<-c("snp.name","pos")

#Get genotypes and SNP names
geno<-geno(vcf)$GT
#snps<-as.data.frame(names(geno))
snps<-as.data.frame(colnames(geno))
rm(vcf)

#Ensure pids in PED file are in the same order as in GENO/TPED
colnames(snps)<-"pids"
ped<-dplyr::left_join(ped,snps, by=c("pid" ="pids"))
tped<-.getTPED(vcf.geno = geno)
map<-.getMAP(tped)
                
########################################################

#5. Save the MAP, TPED, and PED file

filepath.map <- "/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/euro.map"
write.table(map,filepath.map, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)



filepath.tped <- "/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/euro.tped"
write.table(tped,filepath.tped, sep="\t",col.names=FALSE,row.names = TRUE,quote = FALSE)



filepath.ped <- "/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/euro.ped"
write.table(ped,filepath.ped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

########################################################

#6. Run the RV-TDT commands

#OK this works too, what gives?
#Includes './' in the input filename
'/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT' ./NA -G './euro.tped' -P './euro.ped' -M './euro.map' --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

########################################################

# Functions

########################################################

.getTPED<-function(plink = NULL, vcf.geno = NULL){
        
        if(!is.null(plink)){
                tped<-as.data.frame(t(plink[,7:ncol(plink)]))
        } else if (!is.null(vcf.geno)) {
                tped<-as.matrix(splitstackshape::cSplit(vcf.geno, colnames(vcf.geno), c("/")))
                rownames(tped) <-rownames(vcf.geno)
        }
        return(tped)
}

.getMAP<-function(tped){
        
        #TODO: fix the gene.id variable and rename it after the window name
        gene.id<-"NA"
        n.snps<-nrow(tped)
        gene.id.vec<-rep(gene.id,n.snps)
        snps<-rownames(tped)
        mafs<-as.vector(rowMeans(tped))
        map<-as.data.frame(cbind(gene.id.vec,snps,mafs), stringsAsFactors=FALSE)
        return(map)
        
}
