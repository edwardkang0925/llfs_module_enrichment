import os
import re
import ast
from cmath import nan
from typing import List

import pandas as pd
import statsmodels.stats.multitest as smt




def processOnePascalOutput(DIRPATH:str, alpha:float, outputPATH:str):
    """
    Given a path to a pascal output file, extract module index, module genes, and BH-corrected module pvalue

    Args:
        DIRPATH (str): path to a pascal output file
        alpha (float): significance threshold for modules pvalue after BH correction
        outputPATH (str): path to save the whole pascal result

    Returns:
        _type_: list of processed module info, total number of significant pathways
    """
    with open(DIRPATH, "r") as f:
        results = f.read()
        # flatten pval parts
        results = results.replace(",\n", ",")
        results = results.replace(" ", "")
        # 0: module index, 1: module genes, 2: gene uniform pval, 3: module uncorrected pval
        parsedResults = re.findall("\[(.+?),(.*?),array\((.*)\),(.*?)\]", results)
        
        pathwayIndexList = []
        pathwayGenesList = []
        pathwayPvalList = []
        
        for i in range(len(parsedResults)):
            # if a module lost all genes due to missing gene score, exclude it from FDR
            if parsedResults[i][3] != "nan": 
                # ex) "'5'" -> 5
                pathwayIndexList.append(int(parsedResults[i][0].replace("'", "")))
                # string list to list conversion via ast.literal_eval
                pathwayGenesList.append(ast.literal_eval(parsedResults[i][1]))
                pathwayPvalList.append(float(parsedResults[i][3]))
        
        correctedPathwayPvalList = smt.fdrcorrection(pathwayPvalList, alpha)
        
        # output csv file 
        df = pd.DataFrame(list(zip(pathwayIndexList, pathwayGenesList,
                           pathwayPvalList, correctedPathwayPvalList[1])),
                          columns=['moduleIndex', 'moduleGenes', 'modulePval', 'correctedModulePval'])
        df.to_csv(outputPATH)
        
        result = []
        
        numSigPathway = sum(correctedPathwayPvalList[0])
        for tup in zip(pathwayIndexList, pathwayGenesList,
                       correctedPathwayPvalList[0], correctedPathwayPvalList[1]):
            result.append(tup)
        # sort by corrected module pvalue
        result.sort(key=lambda x: x[-1])
        return result, numSigPathway
    
def extractGenesBasedOnPval(DIRPATH:str, pval:float):
    """
    Given a GeneScore file in tsv file format without header and a pvalue threshold, 
    output list of genes with pval less than the threshold

    Args:
        DIRPATH (str): Path to GS file in tsv format. ASSUMPTION: the first col is gene name and the second col is pval
        pval (float): pvalue threshold

    Returns:
        _type_: list of significant genes
    """
    df = pd.read_table(DIRPATH, header=None)
    return list(df[df[1] < pval][0])

    
def recordSignificantModulesFromPascalResult(result, OUTPUTPATH:str, sigGenesList, almostSigGenesList):
    moduleIndexToSize = {}
    moduleIndexSigGenes = {}
    moduleIndexAlmostSigGenes = {}

    for item in result:
        # assumes index of 2 represents bool indicating significance of the module
        if item[2]:
            moduleSizeCounter = 0
            sigGenes = []
            almostSigGenes = []
            with open(OUTPUTPATH.replace(".txt", f"_{item[0]}.txt"), "w") as f:
                for gene in item[1]:
                    if gene in sigGenesList:
                        sigGenes.append(gene)
                    if gene in almostSigGenesList:
                        almostSigGenes.append(gene)
                    f.write(f'{gene}\n')
                    moduleSizeCounter += 1
                    
            moduleIndexToSize[item[0]] = moduleSizeCounter
            moduleIndexSigGenes[item[0]] = sigGenes
            moduleIndexAlmostSigGenes[item[0]] = almostSigGenes
            

    
    return moduleIndexToSize, moduleIndexSigGenes, moduleIndexAlmostSigGenes
                    
                
