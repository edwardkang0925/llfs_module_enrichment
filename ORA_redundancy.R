#!/usr/bin/env Rscript

## TO USE:
## need one folder containing the gene lists of interest
## need another folder containing their reference gene lists
## all gene list files should be single-column .txt
## input-reference list pairs must share the first '.'-separated word of their names
## output files are named based on the first '_'-separated word of the input
## change FIELDS to desired name, method, database, gene id format, input folder, and gene universe folder
## data on the terms with significant overlap will be in SUMMARIES_PATH if such terms exist

library("WebGestaltR")
library(stringr)

setwd("/Users/test/projects/llfs_module_enrichment")

# Things to change
MODULEFILEROOT = "./outputs/parsedPascalOutput/"
SUMMARYROOT = "./outputs/GO_summaries_redundancy"
REPORTROOT = "./outputs/GO_reports_redundancy"

# WebGestalt parameters
METHOD="ORA" # ORA | GSEA | NTA
DATABASE="geneontology_Biological_Process"  # "geneontology_Biological_Process "geneontology_Molecular_Function" # see options with listGeneSet()
GENE_ID="genesymbol" # see options with listIdType()
backgroundSetPATH = './outputs/GOinput/gwas/BMI/gwas_BMI_ppi.txt'
sigmodulePATH = './outputs/parsedPascalOutput/gwas/BMI/significant/gwas_BMI_ppi_2.txt'

setCoverNum = 10
nThreads = 4
filename = paste0("weightedSetCover_", setCoverNum)
apFileName="affinityPropagation_mlogp"


enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                            enrichDatabase=DATABASE, interestGeneFile=sigmodulePATH,
                            interestGeneType=GENE_ID, referenceGeneFile=backgroundSetPATH,
                            referenceGeneType=GENE_ID, isOutput=FALSE, outputDirectory=REPORTROOT)

# WeightedSetCover & Affinity Propagation https://github.com/bzhanglab/WebGestaltR/issues/10
idsInSet <- sapply(enrichResult$overlapId, strsplit, split=";")
names(idsInSet) <- enrichResult$geneSet
minusLogP <- -log(enrichResult$pValue)
minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)

#wscRes <- weightedSetCover(idsInSet, 1 / minusLogP, setCoverNum, nThreads)

apRes <- affinityPropagation(idsInSet, minusLogP)

# make into table (wsc)
weightedGO_full <- enrichResult[c(match(wscRes$topSets, enrichResult$geneSet)),]

# extract representatives from apRes and generate filtered output.
apGO_full <- enrichResult[enrichResult$geneSet %in% apRes$representatives,]

write.csv(weightedGO_full,file.path(SUMMARYROOT,paste0(filename,".csv")),row.names = FALSE)
write.csv(apGO_full,file.path(SUMMARYROOT,paste0(apFileName,".csv")),row.names = FALSE)

