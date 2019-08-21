#Ensembl VEP test

#Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
BiocManager::install("ensemblVEP")

file <- system.file("extdata", "ex2.vcf", package="VariantAnnotation")
file

## When the 'vcf' option is TRUE, a VCF object is returned.

myparam <- VEPFlags(flags=list(vcf=TRUE, host="useastdb.ensembl.org"))
vcf <- ensemblVEP(file, param=myparam)
vcf

## The consequence data are returned as the 'CSQ' column in info.

info(vcf)$CSQ

fp.test.vep <-""
writeVcf()
