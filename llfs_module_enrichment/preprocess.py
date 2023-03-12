import os
import re
import shutil # used for moving files 
from typing import List
import functools as ft


import pandas as pd

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


def queryDirectories(DIRPATH:str) -> List[str]:
    """_summary_
        From a directory, return a list of paths to subdirectories
    
    Args:
        DIRPATH (str): path to the directory where the list of subdirectories are of interest

    Returns:
        List[str]: list of path of subdirectories under DIRPATH 
    """
    result = []
    if os.path.exists(DIRPATH):
        for item in os.listdir(DIRPATH):
            fullPath = os.path.join(DIRPATH, item)
            if os.path.isdir(fullPath):
                result.append(fullPath)
    return result


def querySpecificFiles(DIRPATH:str, endswith='.csv') -> List[str]:
    """
    Given a path to a directory, return list of full path to all csv files under the directory.

    Args:
        DIRPATH (str): Path to the directory containing csv files

    Returns:
        List[str]: List of full path to each of the csv files located under DIRPATH
    """
    result = []
    if os.path.exists(DIRPATH):
        for item in os.listdir(DIRPATH):
            if item.endswith(endswith):
                fullPath = os.path.join(DIRPATH, item)
                result.append(fullPath)
    else:
        print(f"{DIRPATH} does not exist but is queried for CSV files")
    return result

def cmaDirReformat(DIRPATH:str, filename:str, geneNameCol:str, minPvalCol:str) -> None:
    traitDirs = queryDirectories(DIRPATH)
    for traitDir in traitDirs:
        categoryDirs = queryDirectories(traitDir)
        for categoryDir in categoryDirs:
            files = querySpecificFiles(categoryDir)
            for file in files:
                if filename in file:
                    df = pd.read_csv(file)
                    df = df[[geneNameCol, minPvalCol]]
                    csvfilename = os.path.join(traitDir, f"{categoryDir.split('/')[-1]}_{file.split('/')[-1]}")
                    df.to_csv(csvfilename, index=False)
            shutil.rmtree(categoryDir)
            
def gwasDirReformat(DIRPATH:str, geneNameCol:str, pvalCol:str) -> None:
    files = querySpecificFiles(DIRPATH)
    

def combineAcrossCategoriesSelectLowestPval(DIRPATH:str, geneNameCol:str, minPvalCol:str, outputFilePath:str) -> str:
    """
    Given path to a trait folder, merge csv files across categories then retain the minimum pvalue for each gene.

    Args:
        DIRPATH (str): path to trait directory
        geneNameCol (str): gene name column of csv files
        minPvalCol (str): pval column of csv files
        outputFilePath (str): output file path including the filename
    
    Returns:
        str: filepath to the resulting dataframe.
    """
    filesToCombine = querySpecificFiles(DIRPATH)
    dfs = [pd.read_csv(filePATH) for filePATH in filesToCombine]
    # since dfs may have different length depending on the categories, we use outer-merge which takes union of geneNameCol
    # FIXME: suffixes raising some future warning. 
    df_merged = ft.reduce(lambda left,
                          right: pd.merge(left, right, how="outer", 
                                          suffixes=("_x", "_y"), on=geneNameCol), 
                          dfs)
    df_merged[minPvalCol] = df_merged[df_merged.columns[1:]].min(axis=1)
    df_out = df_merged[[geneNameCol, minPvalCol]]
    df_out.drop_duplicates(inplace=True)
    df_out.to_csv(outputFilePath, index=False)
    return outputFilePath
    

def extractCategoryFromFileName(DIRPATH:str, regex:str) -> str:
    """
    Given a path to a csv file, extract the category substring using the given regular expression

    Args:
        DIRPATH (str): path to a csv file containing gene score. 
        regex (str): regular expression to use to extract the category from the file name

    Returns:
        str: _description_
    """
    filename = DIRPATH.split("/")[-1]
    result = re.findall(regex, filename)
    return result[0]

def extractGeneSetFromModuleFile(DIRPATH:str):
    """
    Read a module file and extract a set of genes in the file. 
    It assumes the input file is tsv format where the gene name starts to appear from the thrid column

    Args:
        DIRPATH (str): path to the module file

    Returns:
        _type_: set of genes appear in the module file
    """
    ret = set()
    with open(DIRPATH, "r") as f:
        lines = f.readlines()
        for line in lines:
            columns = line.split()
            for column in columns[2:]:
                ret.add(column)
    return ret
            

def pairwiseProcessGeneScoreAndModule(GSPATH:str, MODULEPATH:str, OUTPUTPATH:str, GOPATH:str, pipeline:str, trait:str, geneNameCol:str, pvalCol:str, sep:str=',') -> None:
    """
    Given a pair of Gene score file and a module file, drop genes which does not exist in either of the files. 
    Write a pair of processed file with the same name. Pascal will take this pair as an input to proceed module enrichment.

    Args:
        GSPATH (str): path to the gene score file
        MODULEPATH (str): path to a pre-defined module file
        OUTPUTPATH (str): path to start building nested subdir for outputs. [pipeline > trait > output]
        pipeline (str): name of the pipeline. e.g. twas, gwas, staar, or cma
        trait (str): name of the trait
        geneNameCol (str): column name for gene name in the GS file
        pvalCol (str): column name for the pvalue in the GS file
        sep (str): if input gene score file is tab separated, pass in '\t' otherwise default should work.
    """
    
    df_gs = pd.read_csv(GSPATH, sep=sep) 
    genesWithScore = set(df_gs[geneNameCol])
    genesInModule = extractGeneSetFromModuleFile(MODULEPATH)
    intersectingGenes = genesWithScore.intersection(genesInModule)
    
    
    # create output directory
    outPvalDir = os.path.join(OUTPUTPATH, pipeline, trait, "pvals")
    createOrCleanDir(outPvalDir, clean=False) # create dir to output processed gs file
    GODir = os.path.join(GOPATH,pipeline, trait)
    createOrCleanDir(GODir, clean=False)
    outModuleDir = os.path.join(OUTPUTPATH, pipeline, trait, "modules")
    createOrCleanDir(outModuleDir, clean=False) # create dir to output processed gs file

    # output files
    gsFileName = GSPATH.split("/")[-1]
    moduleFileName = MODULEPATH.split("/")[-1]
    # output GS file to be used for PASCAL. 
    with open(os.path.join(outPvalDir, f"{pipeline}_{trait}_{moduleFileName[:-4]}.tsv"), "w") as f:
        for index, row in df_gs.iterrows():
            f.write(row[geneNameCol] + "\t" + str(row[pvalCol]) + "\n")
            
    # ouput GO background set 
    with open(os.path.join(GODir, f"{pipeline}_{trait}_{moduleFileName[:-4]}.txt"),"w") as f:
        for index, row in df_gs.iterrows():
            if row[geneNameCol] in intersectingGenes:
                f.write(f"{row[geneNameCol]}\n")
    
    # output module file after intersecting with Gene Score file
    with open(os.path.join(outModuleDir, f"{pipeline}_{trait}_{moduleFileName[:-4]}.tsv"), "w") as f:
        with open(MODULEPATH, "r") as g:
            droppedGeneCounter = 0
            lines = g.readlines()
            for line in lines:
                columns = line.split()
                f.write(columns[0]) # write the module index 
                for column in columns[2:]: # column[1] is always 1.0, so dropped
                    if column in intersectingGenes:
                        f.write("\t" + column)
                    else:
                        droppedGeneCounter += 1
                f.write("\n")
            print(f"Total {droppedGeneCounter} genes were dropped out of {len(genesInModule)} from {moduleFileName}")
            
    
