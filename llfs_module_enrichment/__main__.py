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

# create/clean outputs directories
createOrCleanDir("./outputs/processed/processedGeneScores/")
createOrCleanDir("./outputs/processed/processedModules/")

# get pvals directories (GWAS, TWAS, STAAR or CMA) 
pvalsDirRoot = "./data/pvals/" # location where GWAS, TWAS, STAAR, CMA dirs are

# STAAR with 14 traits and 10 categories.
staar_trait_dirs = queryDirectories(os.path.join(pvalsDirRoot, "staar")) # each trait dir has 10 categories. 
for trait_dir in staar_trait_dirs:
    trait = trait_dir.split("/")[-1]
    print(f"At {trait} directory")
    trait_category_pvals = queryCSVFiles(trait_dir) # each category under a trait 
    for trait_category in trait_category_pvals:
        category = extractCategoryFromFileName(trait_category, f'{trait}_(.*)_staar')  # assumes category substring is between {trait}_ and _staar
        print(f"{category}")
