library(TVTB)

tparam <- TVTBparam(Genotypes(
        ref = "0|0",
        het = c("0|1", "1|0", "0|2", "2|0", "1|2", "2|1"),
        alt = c("1|1", "2|2")),
        ranges = GenomicRanges::GRangesList(
                SLC24A5 = GenomicRanges::GRanges(
                        seqnames = "15",
                        IRanges::IRanges(
                                start = 48413170, end = 48434757)
                )
        )
)

svp <- as(tparam, "ScanVcfParam")
svp

#Step 1: Import phenotypes
phenoFile <- system.file(
        "extdata", "integrated_samples.txt", package = "TVTB")

phenotypes <- S4Vectors::DataFrame(
        read.table(file = phenoFile, header = TRUE, row.names = 1))

phenotypes

#Step 2: Define the VCF file to parse
vcfFile <- system.file(
        "extdata", "chr15.phase3_integrated.vcf.gz", package = "TVTB")
tabixVcf <- Rsamtools::TabixFile(file = vcfFile)

#Step 3: Define VCF import parameters
VariantAnnotation::vcfInfo(svp(tparam)) <- vep(tparam)
VariantAnnotation::vcfGeno(svp(tparam)) <- "GT"
svp(tparam)

#Step 4: Import and pre-process variants

# Import variants as a CollapsedVCF object
vcf <- VariantAnnotation::readVcf(
        tabixVcf, param = tparam, colData = phenotypes)
# Expand into a ExpandedVCF object (bi-allelic records)
vcf <- VariantAnnotation::expand(x = vcf, row.names = TRUE)

################################################################################

initialInfo <- colnames(info(vcf))
vcf <- addOverallFrequencies(vcf = vcf)
setdiff(colnames(info(vcf)), initialInfo)

################################################################################

infoR <- VcfInfoRules(
        exprs = list(
                rare = expression(MAF < 0.01 & MAF > 0),
                common = expression(MAF > 0.05)),
        active = c(TRUE, TRUE)
)
infoR

summary(eval(expr = infoR, envir = vcf))

summary(evalSeparately(expr = infoR, envir = vcf))

vepR <- VcfVepRules(exprs = list(
        missense = expression(Consequence %in% c("missense_variant")),
        CADD_gt15 = expression(CADD_PHRED > 15)
))
vepR
