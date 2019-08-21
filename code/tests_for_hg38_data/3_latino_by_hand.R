#Test 2: Do you still get the error in the Latino dataset if you make the TPED/PED/MAP by hand?

########################################################

setwd("/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets")

########################################################

library(VariantAnnotation)
library(rvtrio)
library(splitstackshape)

########################################################

#1. Read in PED

fp.ped<-"/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/latino.ped"
ped <- read.table(fp.ped, header=FALSE, stringsAsFactors = FALSE)
colnames(ped)<-c("pid","famid","fatid","motid","sex","affected")

########################################################

#2. Read in VCF

filepath.vcf<-"/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/PHASED.8q24.16snp.rarevar.08_19_19.vcf.gz"
hg.assembly<-"hg38"
vcf<-VariantAnnotation::readVcf(filepath.vcf, hg.assembly)

########################################################

#3. Recode VCF

geno(vcf)$GT[geno(vcf)$GT == "0|0"]<-"0/0"
geno(vcf)$GT[geno(vcf)$GT == "0|1"]<-"0/1"
geno(vcf)$GT[geno(vcf)$GT == "1|1"]<-"1/1"
geno(vcf)$GT[geno(vcf)$GT == "1|0"]<-"1/0"

table(geno(vcf)$GT)

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

#Make TPED
map<-.getMAP(tped)

########################################################

#5. Save the MAP, TPED, and PED file

filepath.map <- "/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/latino.map"
write.table(map,filepath.map, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)



filepath.tped <- "/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/latino.tped"
write.table(tped,filepath.tped, sep="\t",col.names=FALSE,row.names = TRUE,quote = FALSE)



filepath.ped <- "/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets/latino.ped"
write.table(ped,filepath.ped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)

########################################################

#6. Run the RV-TD commands

########################################################

#Code below causes seg faults

########################################################

#Without ./ in filename
"/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT" test1 -G latino.map -P latino.ped -M latino.map --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

#With ./ in filename
"/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT" test1 -G ./latino.map -P ./latino.ped -M ./latino.map --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

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
