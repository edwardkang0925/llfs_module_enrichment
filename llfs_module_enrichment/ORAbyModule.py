import os
import pandas as pd
import ast
import subprocess


def outputTxtFile(PATH, filename, genes):
    createOrCleanDir(PATH, False)
    f = open(os.path.join(PATH, filename), 'w')
    for gene in genes:
        f.write(f"{gene}\n")
    f.close()
    
def createOrCleanDir(DIRPATH:str, clean:bool=True) -> None:
    """
    Create or clean a directory

    Args:
        DIRPATH (str): path to the direcotry which needs to be created or cleaned
        clean (bool): flag indicating whether the pre-existing files under the directory should be deleted
    """
    if not os.path.exists(DIRPATH):
        print(f"{DIRPATH} does not exist. Creating ...")
        os.makedirs(DIRPATH)
    else:
        if clean:
            print(f"{DIRPATH} exists. Cleaning up previous files ...")
            for file in os.listdir(DIRPATH):
                os.remove(file)

outputROOTPATH = "/Users/test/projects/llfs_module_enrichment/outputs/goByModules/"
masterSummaryFilePATH = "/Users/test/projects/llfs_module_enrichment/outputs/master_summary.csv"
parsedPascalSummaryROOTPATH = "/Users/test/projects/llfs_module_enrichment/outputs/parsedPascalOutput"
backgroundROOTPATH = "/Users/test/projects/llfs_module_enrichment/outputs/goByModules/background"
sig4genesROOTPATH = "/Users/test/projects/llfs_module_enrichment/outputs/goByModules/sig4genes"

df_summary = pd.read_csv(masterSummaryFilePATH)

df_summary_sig = df_summary[(df_summary.isModuleSig) & (df_summary.network != 'coexpression')]
df_summary_sig = df_summary_sig[(df_summary_sig.study == 'twas') | (df_summary_sig.study == 'cma')]

for index, row in df_summary_sig.iterrows():
    study = row['study']
    trait = row['trait']
    if study == 'twas' and trait == 'mavg_cca':
        trait = 'mavg'
    network = row['network']
    moduleIndex = row["moduleIndex"]
    sig4Genes = ast.literal_eval(row['sig4Genes'])
    pathTORelatedPascalOutput = os.path.join(parsedPascalSummaryROOTPATH, study, trait, f"{study}_{trait}_{network}.csv")
    df = pd.read_csv(pathTORelatedPascalOutput)
    print(pathTORelatedPascalOutput)
    targetModuleGenes = ast.literal_eval(str(df[df.moduleIndex == moduleIndex].moduleGenes.values[0]))
    outputfilename = f"{study}_{trait}_{network}_{moduleIndex}.txt"
    backgroundFileOutPath = os.path.join(backgroundROOTPATH, study, trait)
    sig4genesFileOutPath = os.path.join(sig4genesROOTPATH, study, trait)
    outputTxtFile(backgroundFileOutPath, outputfilename, targetModuleGenes)
    outputTxtFile(sig4genesFileOutPath, outputfilename, sig4Genes)
    
    
subprocess.call("Rscript ./ORAbymodules_redundancy.R", shell=True)