import os
import re
from typing import List
from .preprocess import *
import pandas as pd

def countGOterms(DIRPATH:str)-> int:
    df = pd.read_csv(DIRPATH)
    moduleIndex = DIRPATH.split("/")[-1].split("_")[-1].replace(".csv","")
    return int(moduleIndex), len(df)

def outputMergableORADF(ORAPATH:str, studies:List[str]):
    dict_ora = {'study':[], 'trait':[], 'network':[], 'moduleIndex':[],
                'geneontology_Biological_Process':[],
                'geneontology_Molecular_Function':[]}
    GOtypes = ['geneontology_Biological_Process', 'geneontology_Molecular_Function']
    for study in studies:
        ora_trait_dirs = queryDirectories(os.path.join(ORAPATH, study))
        for ora_trait_dir in ora_trait_dirs:
            trait = ora_trait_dir.split("/")[-1]
            for ora_type in GOtypes:
                module_ora_files = querySpecificFiles(os.path.join(ora_trait_dir, ora_type), endswith=".csv")
                for module_ora_file in module_ora_files:
                    networkType = module_ora_file.split("_")[-2]
                    moduleIndex, GOcount = countGOterms(module_ora_file)
                    # only append network and moduleIndex once for (study, trait) 
                    if ora_type == GOtypes[0]:
                        dict_ora['study'].append(study)
                        dict_ora["trait"].append(trait)
                        dict_ora["network"].append(networkType)
                        dict_ora["moduleIndex"].append(moduleIndex)
                    dict_ora[ora_type].append(GOcount)

    return pd.DataFrame(dict_ora)
