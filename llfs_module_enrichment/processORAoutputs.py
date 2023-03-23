import os
import re
from typing import List
from .preprocess import *
import pandas as pd

def countGOterms(DIRPATH:str)-> int:
    df = pd.read_csv(DIRPATH)
    moduleIndex = DIRPATH.split("/")[-1].split("_")[-1].replace(".csv","")
    maxEnrichmentRatio = -1
    enrichMentRatio = -1
    if len(df) > 0:
        minTermPval = df["FDR"].min()
        maxEnrichmentRatio = df["enrichmentRatio"].max()
        enrichMentRatio = df[df["FDR"] == df["FDR"].min()]["enrichmentRatio"].max()
    else:
        minTermPval = -1
        maxEnrichmentRatio = -1
        enrichMentRatio = -1
    return int(moduleIndex), len(df), minTermPval, enrichMentRatio, maxEnrichmentRatio

def outputMergableORADF(ORAPATH:str, studies:List[str]):
    dict_ora = {'study':[], 'trait':[], 'network':[], 'moduleIndex':[],
                'geneontology_Biological_Process':[], 'BPminCorrectedPval':[],
                'BPminFDREnrichmentRatio':[], 'BPmaxEnrichmentRatio':[],
                'geneontology_Molecular_Function':[], "MFminCorrectedPval":[],
                'MFminFDREnrichmentRatio':[], 'MFmaxEnrichmentRatio':[]}
    GOtypes = ['geneontology_Biological_Process', 'geneontology_Molecular_Function']
    for study in studies:
        ora_trait_dirs = queryDirectories(os.path.join(ORAPATH, study))
        for ora_trait_dir in ora_trait_dirs:
            trait = ora_trait_dir.split("/")[-1]
            for ora_type in GOtypes:
                module_ora_files = querySpecificFiles(os.path.join(ora_trait_dir, ora_type), endswith=".csv")
                for module_ora_file in module_ora_files:
                    networkType = module_ora_file.split("_")[-2]
                    moduleIndex, GOcount, minPval, enrichmentRatio_minFDR, enrichmentRatio_max = countGOterms(module_ora_file)
                    # only append network and moduleIndex once for (study, trait) 
                    if ora_type == GOtypes[0]:
                        dict_ora['study'].append(study)
                        dict_ora["trait"].append(trait)
                        dict_ora["network"].append(networkType)
                        dict_ora["moduleIndex"].append(moduleIndex)
                        dict_ora['BPminCorrectedPval'].append(minPval)
                        dict_ora["BPminFDREnrichmentRatio"].append(enrichmentRatio_minFDR)
                        dict_ora["BPmaxEnrichmentRatio"].append(enrichmentRatio_max)
                    else:
                        dict_ora['MFminCorrectedPval'].append(minPval)
                        dict_ora['MFminFDREnrichmentRatio'].append(enrichmentRatio_minFDR)
                        dict_ora['MFmaxEnrichmentRatio'].append(enrichmentRatio_max)
                    dict_ora[ora_type].append(GOcount)

    return pd.DataFrame(dict_ora)
