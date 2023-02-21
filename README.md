# llfs_module_enrichment

- input
  - pvals
    - GWAS, TWAS, STAAR (with 10 categories), and CMA
    - 14 traits
  - pre-defined modules
    - from Marbach's paper
- output

  - pascal module enrichment result
  - GO analysis file

  preprocess.py contains functions to be used before running Pascal for module enrichment.

postPascal.py contains functions to be used after running Pascal module enrichment
