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
    
    # CMA
    cma_trait_dirs = queryDirectories(os.path.join(pvalsDirRoot, "cma"))
    cmaDirReformat(os.path.join(pvalsDirRoot, "cma"), "CMA_meta", 'markname', 'meta_p')
    for trait_dir in cma_trait_dirs:
        trait = trait_dir.split('/')[-1]
        trait_combined_across_categories = combineAcrossCategoriesSelectLowestPval(trait_dir, geneNameCol="markname", 
                                                                                   minPvalCol="meta_p",
                                                                                   outputFilePath=f"./outputs/log/cma_{trait}_combined.csv")
        for path_to_module_file in os.listdir(PATHTOMODULES):
            if ".txt" in path_to_module_file:
                pairwiseProcessGeneScoreAndModule(trait_combined_across_categories, 
                                                    os.path.join(PATHTOMODULES, path_to_module_file), 
                                                    pathToProcessedInput, GOinputDir,
                                                    "cma", trait, "markname", "meta_p")
             
    # TWAS
    twas_gs_files = querySpecificFiles(os.path.join(pvalsDirRoot, "twas")) # since twas dir has all the csv file, different from staar where each csv files are grouped under a directory <trait> 
    for twas_gs_file in twas_gs_files:
        trait = twas_gs_file.split("_")[1].replace(".csv", "") # HARDCODED location of trait in filename
        for path_to_module_file in os.listdir(PATHTOMODULES):
            if ".txt" in path_to_module_file: # to filter out .DSstore file
                pairwiseProcessGeneScoreAndModule(twas_gs_file, os.path.join(PATHTOMODULES, path_to_module_file),
                                                  pathToProcessedInput, GOinputDir,
                                                  "twas", trait, "Genes", "p_vals_corrected")
    # GWAS
    gwas_gs_files = querySpecificFiles(os.path.join(pvalsDirRoot, 'gwas'))
    for gwas_gs_file in gwas_gs_files:
        trait = gwas_gs_file.split("_")[1].replace(".csv","") # HARDCODED
        for path_to_module_file in os.listdir(PATHTOMODULES):
            if ".txt" in path_to_module_file: # to filter out .DSstore file 
                pairwiseProcessGeneScoreAndModule(gwas_gs_file, os.path.join(PATHTOMODULES, path_to_module_file),
                                                  pathToProcessedInput, GOinputDir,
                                                  "gwas", trait, "Genes", "p_vals")
    
else:
    geneScoreDir = "./outputs/pascalInput/"
    pascalOutputDir = "./outputs/pascalOutput/"
    OUTPUTDIR = "./outputs/parsedPascalOutput/"
    MASTER_SUMMARY_OUTPATH = "./outputs/master_summary.csv"
    ORAPATH = "./outputs/GO_summaries/"
    ora_types = ['geneontology_Biological_Process', 'geneontology_Molecular_Function']
    ORA_SUMMARY_PATH = "./outputs/ora_summary.csv"
    studies = ['staar', 'twas', 'gwas', 'cma'] # dir name
    NUMTWASGENES = 17972
    NUMGWASGENES = 23534
    sigPvalThreshold = {
        "staar-adjTC": 0.05 / (145529 - 10),
        "staar-fev1fvc": 0.05 / (145057 - 10),
        "staar-adjLDLF": 0.05 / (145503 - 10),
        "staar-BMI": 0.05 / (146004 - 10),
        "staar-pulse": 0.05 / (145971 - 10),
        "staar-fhshdl": 0.05 / (146056 - 10),
        "staar-lnTG": 0.05 / (145583 - 10),
        "staar-ABI": 0.05 / (142356 - 10),
        "staar-waist": 0.05 / (145890 - 10),
        "staar-fvc": 0.05 / (145061 - 10),
        "staar-fev1": 0.05 / (145122 - 10),
        "cma-adjTC": 0.05 / (177926 - 10),
        "cma-fev1fvc": 0.05 / (177770 - 10),
        "cma-adjLDLF": 0.05 / (177923 - 10),
        "cma-BMI": 0.05 / (178072 - 10),
        "cma-pulse": 0.05 / (178059 - 10),
        "cma-fhshdl": 0.05 / (178116 - 10),
        "cma-lnTG": 0.05 / (177943 - 10),
        "cma-ABI": 0.05 / (176861 - 10),
        "cma-waist": 0.05 / (178037 - 10),
        "cma-fvc": 0.05 / (177758 - 10),
        "cma-fev1": 0.05 / (177788 - 10),
        "gwas": 0.05 / NUMGWASGENES,
        "twas": 0.05 / NUMTWASGENES
    }
        
    # master summary file columns
    summary_dict = {'study':[],
                    'trait':[],
                    'network':[],
                    'moduleIndex':[],
                    "isModuleSig":[],
                    "modulePval":[],
                    "moduleBonPval":[],
                    'size':[],
                    'numSigGenes':[],
                    'sigGenes':[],
                    'sig1Genes':[],
                    'sig2Genes':[],
                    'sig3Genes':[],
                    'sig4Genes':[]
                    }
    
    for study in studies:
        pascal_trait_dirs = queryDirectories(os.path.join(pascalOutputDir, study))
        for pascal_trait_dir in pascal_trait_dirs:
            pascalOutputFiles = querySpecificFiles(pascal_trait_dir, endswith='.txt')
            trait = pascal_trait_dir.split("/")[-1]
            if study == "cma" or study == "staar":
                studyCode = f"{study}-{trait}"
            else:
                studyCode = study
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
                                                                     pascalOutputFileName.replace(".txt", ".tsv")), sigPvalThreshold[studyCode])
                sig1GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', 
                                                                     pascalOutputFileName.replace(".txt", ".tsv")), sigPvalThreshold[studyCode]*10)
                sig2GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', 
                                                                     pascalOutputFileName.replace(".txt", ".tsv")), sigPvalThreshold[studyCode]*100)
                sig3GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', 
                                                                     pascalOutputFileName.replace(".txt", ".tsv")), sigPvalThreshold[studyCode]*1000)
                sig4GenesList = extractGenesBasedOnPval(os.path.join(geneScoreDir, study, trait, 'pvals', 
                                                                     pascalOutputFileName.replace(".txt", ".tsv")), sigPvalThreshold[studyCode]*10000)
                moduleToSize, moduleToPval, moduleToCorrectedPval, isModuleSig, sigGenesDict, sig1GenesDict, sig2GenesDict, sig3GenesDict, sig4GenesDict = recordModulesFromPascalResult(result, os.path.join(sigModuleOutPath, pascalOutputFileName),
                                                                                                          sigGenesList, sig1GenesList, sig2GenesList, sig3GenesList, sig4GenesList) 
                
                # Master summary file data
                for moduleIndex in sigGenesDict.keys():
                    summary_dict['study'].append(study)
                    summary_dict['trait'].append(trait)
                    summary_dict['network'].append(networkType)
                    summary_dict['moduleIndex'].append(moduleIndex)
                    summary_dict['isModuleSig'].append(isModuleSig[moduleIndex])
                    summary_dict['modulePval'].append(moduleToPval[moduleIndex])
                    summary_dict['moduleBonPval'].append(moduleToCorrectedPval[moduleIndex])
                    summary_dict['size'].append(moduleToSize[moduleIndex])
                    summary_dict['numSigGenes'].append(len(sigGenesDict[moduleIndex]))
                    summary_dict['sigGenes'].append(sigGenesDict[moduleIndex])
                    summary_dict["sig1Genes"].append(sig1GenesDict[moduleIndex])
                    summary_dict['sig2Genes'].append(sig2GenesDict[moduleIndex])
                    summary_dict['sig3Genes'].append(sig3GenesDict[moduleIndex])
                    summary_dict['sig4Genes'].append(sig4GenesDict[moduleIndex])

    # output summary file
    df_summary = pd.DataFrame(summary_dict)
    
    # Run GO enrichment and output summary
    subprocess.call("Rscript ./ORA_redundancy_batch.R", shell=True)
    df_ora = outputMergableORADF(ORAPATH, studies)
    
    # HARDCODED for traitname with underscore
    df_ora['trait'] = df_ora['trait'].replace("mavg", "mavg_cca")
    df_summary['trait'] = df_summary['trait'].replace('mavg', 'mavg_cca')
    df_merge = pd.merge(df_summary, df_ora, how='left', on=['study','trait','network', 'moduleIndex'])
    df_merge.fillna(-1, inplace=True)
    df_merge = df_merge.astype({"geneontology_Biological_Process":'int', "geneontology_Molecular_Function": 'int'})
    df_merge.sort_values(by=['isModuleSig', 'numSigGenes'], inplace=True, ascending=False)
    df_merge.to_csv(MASTER_SUMMARY_OUTPATH)