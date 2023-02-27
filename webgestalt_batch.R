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

setwd("/Users/test/projects/llfs_module_enrichment")

# PROJECT_NAME names summary folder, reports folders will be "Project_"NAME_METHOD
PROJECT_NAME="test_run_02272023"
METHOD="ORA" # ORA | GSEA | NTA
DATABASE="geneontology_Biological_Process" # "geneontology_Biological_Process "geneontology_Molecular_Function" # see options with listGeneSet()
GENE_ID="genesymbol" # see options with listIdType()

INPUT_PATH="./outputs/parsedPascalOutput/staar/ABI/significant" # path to folder of module genes
REFERENCE_PATH="./outputs/GOinput/staar/ABI" # path to folder of background gene lists

# reports are more in-depth than summaries - advisable to keep reports FALSE if not needed
REPORTS_PATH="./outputs/GO_results" # only used if GENERATE_REPORT=TRUE
SUMMARIES_PATH="./outputs/GO_summaries" # will be created if does not exist
GENERATE_REPORT=FALSE

# path must exist even if GENERATE_REPORT=FALSE
if (!dir.exists(REPORTS_PATH)) {
  dir.create(REPORTS_PATH, recursive=TRUE)
}

if (!dir.exists(file.path(SUMMARIES_PATH, PROJECT_NAME))) {
  dir.create(file.path(SUMMARIES_PATH, PROJECT_NAME), recursive=TRUE)
}


if (dir.exists(INPUT_PATH) & dir.exists(REFERENCE_PATH)) {
  # used to pair input with reference
  ref_list = list.files(REFERENCE_PATH)
  ref_list_idents = lapply(ref_list, function(x) unlist(strsplit(x, '[.]'))[1]) # get list of module index
  
  # used to allow for duplicate names
  prev_done = c()
  
  for (file_name in list.files(INPUT_PATH)) {
    # get name for input file
    name = unlist(strsplit(file_name, '.txt'))[1]
    
    tf_method = paste0(name, '_', METHOD)
    
    prev_done = c(prev_done, name)
    
    # find proper reference list
    ident = unlist(strsplit(file_name, '[.]'))[1]
    ref_index = which(ident == ref_list_idents)[1]
    
    # perform enrichment analysis
    enrich_df <- WebGestaltR(
      enrichMethod = METHOD,
      organism = "hsapiens",
      enrichDatabase = DATABASE,
      interestGeneFile = file.path(INPUT_PATH, file_name),
      interestGeneType = GENE_ID,
      referenceGeneFile = file.path(REFERENCE_PATH, list.files(REFERENCE_PATH)[ref_index]),
      referenceGeneType = GENE_ID,
      minNum = 10, # default 10
      maxNum = 500, # default 500
      reportNum = 20, # default 20
      isOutput = GENERATE_REPORT,
      outputDirectory = REPORTS_PATH,
      projectName = tf_method
    )
    
    # save summary as a .csv file
    if (!is.null(enrich_df)) {
      sig_df <- subset(enrich_df, select = -c(link))
      if (nrow(sig_df) > 0) {
        sig_df['database'] <- rep(DATABASE, nrow(sig_df))
        write.csv(sig_df,file.path(SUMMARIES_PATH,PROJECT_NAME,paste0(tf_method,"_summary.csv")),row.names = FALSE)
      } else {
        print("NO SIGNIFICANT OVERLAPS")
        write.csv(NULL,file.path(SUMMARIES_PATH,PROJECT_NAME,paste0(tf_method,"_summary.csv")),row.names = FALSE)
      }
    } else {
      print("NO SIGNIFICANT OVERLAPS")
      write.csv(NULL,file.path(SUMMARIES_PATH,PROJECT_NAME,paste0(tf_method,"_summary.csv")),row.names = FALSE)
      
    }
  }
} else {
  print("INPUT NOT FOUND")
}
