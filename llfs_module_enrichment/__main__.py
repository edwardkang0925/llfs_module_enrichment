import sys
import os
from importlib.metadata import version

from .preprocess import *

'''
Before running this pipeline,
TODO
    * place all the pvals directories under data/pvals
'''


def main():
    print(version('llfs_module_enrichment'))

# get pvals directories (GWAS, TWAS, STAAR or CMA) 
pvalsDirRoot = "./data/pvals/" # location where GWAS, TWAS, STAAR, CMA dirs are
PATHTOMODULES = "./data/modules/cherryPickModules/"
pathToProcessedInput = "./outputs/pascalInput/"

# STAAR with 14 traits and 10 categories.
staar_trait_dirs = queryDirectories(os.path.join(pvalsDirRoot, "staar")) # each trait dir has 10 categories. 
for trait_dir in staar_trait_dirs:
    trait = trait_dir.split("/")[-1]
    print(f"At {trait} directory")
    trait_category_pval_files = queryCSVFiles(trait_dir) # each category under a trait 
    for trait_category_file in trait_category_pval_files:
        category = extractCategoryFromFileName(trait_category_file, f'{trait}_(.*)_staar')  # assumes category substring is between {trait}_ and _staar
        print(f"{category}")
        for path_to_module_file in os.listdir(PATHTOMODULES):
            if ".txt" in path_to_module_file:
                pairwiseProcessGeneScoreAndModule(trait_category_file, os.path.join(PATHTOMODULES, path_to_module_file), pathToProcessedInput, 
                                                  "staar", trait, "hgnc_symbol", "pval")
            
# TWAS with 11 traits without category
twas_gs_files = queryCSVFiles(os.path.join(pvalsDirRoot, "twas")) # since twas dir has all the csv file, different from staar where each csv files are grouped under a directory <trait> 
for twas_gs_file in twas_gs_files:
    trait = twas_gs_file.split("_")[7] # HARDCODED location of trait in filename
    for path_to_module_file in os.listdir(PATHTOMODULES):
        if ".txt" in path_to_module_file:
            pairwiseProcessGeneScoreAndModule(twas_gs_file, os.path.join(PATHTOMODULES, path_to_module_file), pathToProcessedInput, 
                                                "twas", trait, "HGNC", "p_vals_corrected")
