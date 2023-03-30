library(WebGestaltR)
library(biomaRt)

setwd("/Users/test/projects/llfs_module_enrichment")


listGeneSet(organism="hsapiens")

ensembl=useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')
getBM(attributes = c('hgnc_symbol'), 
      filters = 'go', 
      values = 'GO:0030098', 
      mart = ensembl)


bp_database <- loadGeneSet(organism='hsapiens', enrichDatabase = "geneontology_Biological_Process")
geneSet <- bp_database$geneSet # each row is a go term. gene column has int gene id. Seems I need to use mapUserId() to get member genes. 

# [loadInterestGene](https://github.com/bzhanglab/WebGestaltR/blob/55d47be3843735cc9bf90f76706026d8a68a98e0/R/loadGeneList.R)

# [mapUserId](https://github.com/bzhanglab/WebGestaltR/blob/5c75381918ce9e2bf78f9d4a0f9d66d118ad02b8/R/reportUtils.R)

GENE_ID="genesymbol" # see options with listIdType()

map <- idMapping(organism='hsapiens',
                dataType='list', 
                inputGeneFile='./outputs/log/allGenes.txt',
                sourceIdType=GENE_ID, 
                targetIdType=NULL, 
                mappingOutput=FALSE)

# output list of unmapped genes. 
fileConn <- file("unmappedGenes.txt")
writeLines(map$unmapped, fileConn)
close(fileConn)

