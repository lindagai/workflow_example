#Test 1: Can you reproduce the error in the european datasets?

########################################################

setwd("/Users/lindagai 1/Documents/classes/4th year/Research/latino_datasets")

########################################################

library(VariantAnnotation)
library(rvtrio)
library(splitstackshape)

########################################################

#1. Read in PED and VCF

filepath.vcf<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/data/8q24.cleaned.phased.filtered.annotation.rarevar.monomorphs.removed.recode.vcf"
filepath.vcf.ped<-"/Users/lindagai 1/Documents/classes/4th year/Research/rvtrio/data/gmkf_euro_completetrios_ids_match_vcf_mend_errors_removed.phen"

ped<-read.table(filepath.vcf.ped,header=FALSE)
colnames(ped)<-c("pid","famid","fatid","motid","sex","affected")
head(ped)

vcf<-VariantAnnotation::readVcf(filepath.vcf, "hg19")

########################################################

#2. Run RV_TDT and keep output files instead of deleting them

RV_TDT.results<-RV_TDT(vcf=vcf, vcf.ped = ped, rv.tdt.dir = filepath.to.RV_TDT)
RV_TDT.results

########################################################

#3. Run RV_TDT on output files with command made by wrapper functions

#This is the command that works
#Includes './' in the input filename
'/Users/lindagai 1/Documents/classes/4th year/Research/rv-tdt-master/rvTDT' ./NA -G './window1.1-368M.tped' -P './pedfile.ped' -M './window1.1-368M.map' --adapt 500 --alpha 1e-05 --permut 2000 --lower_cutoff 0 --upper_cutoff 100 --minVariants 3 --maxMissRatio 1

########################################################

# RV_TDT edited to print RV-TDT commands

########################################################

RV_TDT<-function(plink.ped=NULL, vcf = NULL, vcf.ped = NULL, rv.tdt.dir, window.size=0, window.type = "M", adapt = 500, alpha = 0.00001, permut = 2000, lower_cutoff = 0, upper_cutoff = 100, minVariants = 3, maxMissRatio = 1){
        
        #Extract parameters
        parameters<-c(adapt, alpha, permut,lower_cutoff, upper_cutoff, minVariants,maxMissRatio)
        
        #Obtain ped/tped and df of snp.names/positions based on file input type
        if(!is.null(plink.ped)){
                ped <-plink.ped[,1:6]
                tped<-.getTPED(plink = plink.ped)
                
                #Ensure pids in PED file are in the same order as in GENO/TPED
        } else if (!is.null(vcf) & !is.null(vcf.ped)){
                #Get positions and SNP names
                snp.pos.df<-as.data.frame(cbind(names(vcf),start(rowRanges(vcf))))
                colnames(snp.pos.df)<-c("snp.name","pos")
                
                #Get genotypes and SNP names
                #TODO: Is the as.data.frame necessary?
                geno<-geno(vcf)$GT
                #snps<-as.data.frame(names(geno))
                snps<-as.data.frame(colnames(geno))
                rm(vcf)
                
                #Ensure pids in PED file are in the same order as in GENO/TPED
                colnames(snps)<-"pids"
                ped<-dplyr::left_join(ped,snps, by=c("pid" ="pids"))
                tped<-.getTPED(vcf.geno = geno)
        } else {
                print("Check your input files.")
                return(NULL)
        }
        
        #Get input files
        map<-.getMAP(tped)
        
        results<-.runRV_TDT(ped, map, tped, rv.tdt.dir, window.size, snp.pos.df, param = parameters)
        return(results)
        
}

########################################################

# Helper Functions - RV_TDT

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

############################

#1. gene
#2. variant id. The variant id must matches with the variant id in tped file
#3. maf

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

############################

.runRV_TDT<-function(ped, map, tped, rv.tdt.dir, window.size=0, snp.pos.df, param, window.type  = "M"){
        
        n.snps<- nrow(map)
        
        #return rv_tdt results as df
        if (window.size==0){
                n.windows<-1
                window.size<-n.snps
        } else {
                n.windows<- max(n.snps - window.size + 1,1)
        }
        
        results.df<-as.data.frame(matrix(data=NA,nrow=n.windows, ncol=10))
        colnames(results.df)<-c(
                "gene.name",
                "CMC.Analytical","BRV.Haplo","CMC.Haplo","VT.BRV.Haplo","VT.CMC.Haplo","WSS.Haplo",
                "start.pos","mid.window.pos","end.pos"
        )
        
        #TODO: Get rid of for loop
        #Make this into a function and use apply instead
        #Indexing is very slow for data frames!!
        
        for (i in (1:n.windows)){
                start.index<-i
                end.index<-min(i+window.size-1,n.snps)
                mid.index<- floor(start.index + (end.index - start.index)/2)
                input.filepaths<-.getInputFilesForWindow(ped,map,tped,window.type,start.index,end.index,i)
                curr.window.result<-.runRV_TDTOnWindow(input.filepaths,rv.tdt.dir,param)
                pos.info<-.getWindowPos(start.index,mid.index,end.index,snp.pos.df)
                results.df[i,1:7]<-curr.window.result
                results.df[i,8:10]<-pos.info
                
        }
        
        return(results.df)
        
}

########################################################

# Helper Functions - runRV_TDT

########################################################

.getInputFilesForWindow<-function(ped, map, tped, window.type,start.index,end.index,i){
        
        #TODO: fix this to get the directory functions.R is in
        data.dir<-"./"
        
        file.param<-paste0("window",i,".",
                           start.index,"-",end.index,window.type)
        file.param
        
        filepath.map.new<-paste0(data.dir,file.param,".map")
        filepath.tped.new<-paste0(data.dir,file.param,".tped")
        filepath.tped.new
        filepath.ped<-paste0(data.dir,"pedfile.ped")
        filepath.ped
        
        sm.map<-map[start.index:end.index,]
        sm.tped<-tped[start.index:end.index,]
        
        write.table(sm.map,filepath.map.new, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        write.table(sm.tped,filepath.tped.new, sep="\t",col.names=FALSE,row.names = TRUE,quote = FALSE)
        
        if (i==1){
                #Make sure ped file sample IDs are the same order as tped
                head(ped)
                write.table(ped,filepath.ped, sep="\t",col.names=FALSE,row.names = FALSE,quote = FALSE)
        }
        
        filepaths<-c(filepath.tped.new, filepath.ped, filepath.map.new)
}

############################

.runRV_TDTOnWindow<-function(input.filepaths,rv.tdt.dir,param){
        
        .calculateRV_TDTOnWindow(input.filepaths,rv.tdt.dir, param)
        results<-.extract.results()
        #.clean.up.rv_tdt(input.filepaths)
        return(results)
        
}

############################

.calculateRV_TDTOnWindow<-function(input.filepaths,rv.tdt.dir,param){
        
        #TODO: Fix u to be 0.01
        #TODO: Add in all the parameters that RV-TDT includes
        #TODO: Modify command accordingly
        
        adapt<-param[1]
        
        alpha<-param[2]
        
        permut<-param[3]
        
        lower_cutoff<-param[4]
        
        upper_cutoff<-param[5]
        
        minVariants<-param[6]
        
        maxMissRatio <-param[7]
        
        filepath.tped<-paste0("'",input.filepaths[1],"'")
        filepath.phen<-paste0("'",input.filepaths[2],"'")
        filepath.map<-paste0("'",input.filepaths[3],"'")
        rv.tdt.dir<-paste0("'",rv.tdt.dir,"'")
        
        #filepath.tped<-paste0(input.filepaths[1])
        #filepath.phen<-paste0(input.filepaths[2])
        #filepath.map<-paste0(input.filepaths[3])
        
        gene.name<-"NA"
        
        #TODO: Add directories for each window?
        rv.tdt.results.dir<-paste0("./",gene.name)
        rv.tdt.results.dir
        
        command<-paste0(rv.tdt.dir, " ", rv.tdt.results.dir,
                        " -G ", filepath.tped,
                        " -P ", filepath.phen,
                        " -M ", filepath.map,
                        " --adapt ", adapt,
                        " --alpha ", alpha,
                        " --permut ", permut,
                        " --lower_cutoff ", lower_cutoff,
                        " --upper_cutoff ", upper_cutoff,
                        " --minVariants ", minVariants,
                        " --maxMissRatio ", maxMissRatio
        )
        
        print(command)
        system(command)
        
}

############################

.getWindowPos<-function(start.index,mid.index,end.index,snp.pos.df){
        
        # Get start, end, middle position of each window
        start.pos<-as.character(snp.pos.df$pos[start.index])
        mid.pos<-as.character(snp.pos.df$pos[start.index])
        end.pos<-as.character(snp.pos.df$pos[end.index])
        pos.info<-c(start.pos,mid.pos,end.pos)
        
        # Add to results df for easier graphing
        return(pos.info)
        
}

########################################################

# Helper Functions - .calculateRV_TDTOnWindow

########################################################

.extract.results<-function(){
        gene.name<-"NA"
        filepath.results<-paste0("./",gene.name,"_pval/",gene.name,".pval")
        filepath.results
        pval.df<-read.table(filepath.results,comment.char = "",header=TRUE)
        pval.df
        return(pval.df)
        
}

############################

.clean.up.rv_tdt<-function(input.filepaths){
        
        #delete all input files
        filepath.tped<-paste0("'",input.filepaths[1],"'")
        filepath.map<-paste0("'",input.filepaths[3],"'")
        
        delete.input.files<-paste("rm",
                                  filepath.tped,
                                  filepath.map
        )
        
        #print(delete.input.files)
        system(delete.input.files)
        
        #delete all output files
        currwd<-getwd()
        filepath.pval<-paste0("'",currwd,"/NA_pval/NA.pval","'")
        filepath.results<-paste0("'",currwd,"/NA_rvTDT/NA.rvTDT","'")
        delete.output.files<-paste("rm", filepath.pval, filepath.results)
        #print(delete.output.files)
        system(delete.output.files)
        
        #delete all output directories
        dir.pval<-paste0("'",currwd,"/NA_pval","'")
        dir.results<-paste0("'",currwd,"/NA_rvTDT","'")
        delete.output.dir<-paste("rmdir", dir.pval, dir.results)
        #print(delete.output.dir)
        system(delete.output.dir)
        
}


