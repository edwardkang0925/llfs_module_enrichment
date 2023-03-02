import sys
import os
import argparse
import subprocess
from importlib.metadata import version

from .preprocess import *
from .postPascal import *
from .processORAoutputs import *
'''
Before running this pipeline,
TODO
    * place all the pvals directories under data/pvals
    * outputs dir better be emptied especially the pascalInput dir
'''

def main():
    print(version('llfs_module_enrichment'))


PREPROCESS = False

if PREPROCESS:
    # get pvals directories (GWAS, TWAS, STAAR or CMA) 
    pvalsDirRoot = "./data/pvals/" # location where GWAS, TWAS, STAAR, CMA dirs are
    PATHTOMODULES = "./data/modules/cherryPickModules/"
    pathToProcessedInput = "./outputs/pascalInput/"
    GOinputDir = "./outputs/GOinput/"
    createOrCleanDir('./outputs/log/') # where combined staar file will be saved
    
    # STAAR with 14 traits and take min pval across categories
    staar_trait_dirs = queryDirectories(os.path.join(pvalsDirRoot, "staar")) # each trait dir has 10 categories. 
    for trait_dir in staar_trait_dirs:
        trait = trait_dir.split("/")[-1]
        print(f"At {trait} directory")
        trait_combined_across_categories = combineAcrossCategoriesSelectLowestPval(trait_dir, geneNameCol="hgnc_symbol", 
                                                                                   minPvalCol="pval",
                                                                                   outputFilePath=f"./outputs/log/{trait}_combined_staar_s.csv")
        for path_to_module_file in os.listdir(PATHTOMODULES):
            if ".txt" in path_to_module_file: # to filter out .DSstore file 
                pairwiseProcessGeneScoreAndModule(trait_combined_across_categories, 
                                                    os.path.join(PATHTOMODULES, path_to_module_file), 
                                                    pathToProcessedInput, GOinputDir,
                                                    "staar", trait, "hgnc_symbol", "pval")
    
                
    # TWAS with 11 traits without category
    twas_gs_files = querySpecificFiles(os.path.join(pvalsDirRoot, "twas")) # since twas dir has all the csv file, different from staar where each csv files are grouped under a directory <trait> 
    for twas_gs_file in twas_gs_files:
        trait = twas_gs_file.split("_")[7] # HARDCODED location of trait in filename
        for path_to_module_file in os.listdir(PATHTOMODULES):
            if ".txt" in path_to_module_file: # to filter out .DSstore file 
                pairwiseProcessGeneScoreAndModule(twas_gs_file, os.path.join(PATHTOMODULES, path_to_module_file),
                                                  pathToProcessedInput, GOinputDir,
                                                  "twas", trait, "HGNC", "p_vals_corrected")
else:
    geneScoreDir = "./outputs/pascalInput/"
    pascalOutputDir = "./outputs/pascalOutput/"
    OUTPUTDIR = "./outputs/parsedPascalOutput/"
    MASTER_SUMMARY_OUTPATH = "./outputs/master_summary.csv"
    ORAPATH = "./outputs/GO_summaries/"
    ora_types = ['geneontology_Biological_Process', 'geneontology_Molecular_Function']
    ORA_SUMMARY_PATH = "./outputs/ora_summary.csv"
    studies = ['staar', 'twas'] # dir name
    sigPvalThreshold = {'staar':2.5*(10**-7), 'twas':2.5*(10**-6)}
    almostSigPvalThreshold = {'staar':2.5*(10**-5), 'twas':2.5*(10**-5)}
    
    # master summary file columns
    summary_dict = {'study':[],
                    'trait':[],
                    'network':[],
                    'moduleIndex':[],
                    'size':[],
                    'numSigGenes':[],
                    'sigGenes':[],
                    'numAlmostSigGenes':[],
                    'almostSigGenes':[]
                    }
    
    for study in studies:
        pascal_trait_dirs = queryDirectories(os.path.join(pascalOutputDir, study))
        for pascal_trait_dir in pascal_trait_dirs:
            pascalOutputFiles = querySpecificFiles(pascal_trait_dir, endswith='.txt')
            trait = pascal_trait_dir.split("/")[-1]
            for pascalOutputFile in pascalOutputFiles:
                pascalOutputFileName = pascalOutputFile.split("/")[-1]
                networkType = pascalOutputFile.split("_")[-1].replace(".txt", "")
                outputpath = os.path.join(OUTPUTDIR,study,trait)
                createOrCleanDir(outputpath, clean=False)
                result, numSigModules = processOnePascalOutput(pascalOutputFile, alpha=0.05,
                                                               outputPATH=os.path.join(outputpath, pascalOutputFileName.replace(".txt",".csv")))  
                sigModuleOutPath = os.path.join(outputpath, "significant")
                createOrCleanDir(sigModuleOutPath, clean=False)
                sigGenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', 
                                                                     pascalOutputFileName.replace(".txt", ".tsv")), sigPvalThreshold[study])
                almostSigGenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', 
                                                                     pascalOutputFileName.replace(".txt", ".tsv")), almostSigPvalThreshold[study])

                moduleToSize, sigGenesDict, almostSigGenesDict = recordSignificantModulesFromPascalResult(result, os.path.join(sigModuleOutPath, pascalOutputFileName),
                                                                                                          sigGenesList, almostSigGenesList) 
                
                # Master summary file data
                for moduleIndex in sigGenesDict.keys():
                    summary_dict['study'].append(study)
                    summary_dict['trait'].append(trait)
                    summary_dict['network'].append(networkType)
                    summary_dict['moduleIndex'].append(moduleIndex)
                    summary_dict['size'].append(moduleToSize[moduleIndex])
                    summary_dict['numSigGenes'].append(len(sigGenesDict[moduleIndex]))
                    summary_dict['sigGenes'].append(sigGenesDict[moduleIndex])
                    summary_dict['numAlmostSigGenes'].append(len(almostSigGenesDict[moduleIndex]))
                    summary_dict["almostSigGenes"].append(almostSigGenesDict[moduleIndex])
                                                             
    
    # output summary file
    df_summary = pd.DataFrame(summary_dict)
    # HARDCODED for traitname with underscore
    df_summary['trait'] = df_summary['trait'].replace('mavg', 'mavg_cca')
    
    
    # Run GO enrichment and output summary
    # FIXME: after making ora_summary dataframe, instead of outputting it, merge to master summary file.
    #subprocess.call("Rscript ./webgestalt_batch.R", shell=True)
    
    df_ora = outputMergableORADF(ORAPATH, studies)
    # HARDCODED for traitname with underscore
    df_ora['trait'] = df_ora['trait'].replace("mavg", "mavg_cca")
    df_merge = pd.merge(df_summary, df_ora, how='left', on=['study','trait','network', 'moduleIndex'])
    df_merge.to_csv(MASTER_SUMMARY_OUTPATH)
    

            
    
