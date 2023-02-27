import sys
import os
import argparse
from importlib.metadata import version

from .preprocess import *
from .postPascal import *
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
    SUMMARYOUTPATH = "./outputs/master_summary.csv"
    studies = ['staar', 'twas'] # dir name
    sigPvalThreshold = {'staar':2.5*(10**-7), 'twas':2.5*(10**-6)}
    almostSigPvalThreshold = {'staar':2.5*(10**-5), 'twas':2.5*(10**-5)}
    
    # master summary file columns
    studyList = []
    traitList = []
    networkList = []
    numSigModuleList = []
    moduleIndexToModuleSizeList = [] # will be list of list. ex) [[10,20,30], [11,13]]
    moduleIndexToSig7GenesList = []
    moduleIndexToSig6GenesList = []
    moduleIndexToSig5GenesList = []
    moduleIndexToSig4GenesList = []
    moduleIndexToSig3GenesList = []
    moduleIndexToSig2GenesList = []
    numSigGenesInSigModulesList = []
    
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
                print(f"{pascalOutputFileName} has {numSigModules} enriched modules")
                sigModuleOutPath = os.path.join(outputpath, "significant")
                createOrCleanDir(sigModuleOutPath, clean=False)
                sig7GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', pascalOutputFileName.replace(".txt", ".tsv")), 2.5*(10**-7))
                sig6GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', pascalOutputFileName.replace(".txt", ".tsv")), 2.5*(10**-6))
                sig5GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', pascalOutputFileName.replace(".txt", ".tsv")), 2.5*(10**-5))
                sig4GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', pascalOutputFileName.replace(".txt", ".tsv")), 2.5*(10**-4))
                sig3GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', pascalOutputFileName.replace(".txt", ".tsv")), 2.5*(10**-3))
                sig2GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', pascalOutputFileName.replace(".txt", ".tsv")), 2.5*(10**-2))


                moduleToSize, sig7, sig6, sig5, sig4, sig3, sig2 = recordSignificantModulesFromPascalResult(result, os.path.join(sigModuleOutPath, pascalOutputFileName),
                                                                                                        sig7GenesList, sig6GenesList, sig5GenesList, sig4GenesList, sig3GenesList, sig2GenesList)
                
                # Master summary file dat
                studyList.append(study)
                traitList.append(trait)
                networkList.append(networkType)
                numSigModuleList.append(numSigModules)
                moduleIndexToModuleSizeList.append(moduleToSize)
                moduleIndexToSig7GenesList.append(sig7)
                moduleIndexToSig6GenesList.append(sig6)
                moduleIndexToSig5GenesList.append(sig5)
                moduleIndexToSig4GenesList.append(sig4)
                moduleIndexToSig3GenesList.append(sig3)
                moduleIndexToSig2GenesList.append(sig2)
                if study == "staar":
                    numSigGenesInSigModulesList.append(sum([len(l) for l in sig7.values()]))
                else:
                    numSigGenesInSigModulesList.append(sum([len(l) for l in sig2.values()]))

        df_summary = pd.DataFrame(list(zip(studyList, traitList, networkList, numSigModuleList, moduleIndexToModuleSizeList,
                                   moduleIndexToSig7GenesList, moduleIndexToSig6GenesList, moduleIndexToSig5GenesList,
                                   moduleIndexToSig4GenesList, moduleIndexToSig3GenesList, moduleIndexToSig2GenesList,
                                   numSigGenesInSigModulesList)),
                                  columns=['study', 'trait', 'network', 'numSigModules', 'moduleIndexToSize', 'moduleIndexToSigGenes7', 'moduleIndexToSigGenes6',
                                           'moduleIndexToSigGenes5','moduleIndexToSigGenes4', 'moduleIndexToSigGenes3', 'moduleIndexToSigGenes2', 'numSigGenesInSigModules'])
        df_summary.to_csv(SUMMARYOUTPATH)