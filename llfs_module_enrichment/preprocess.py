import os
import re
from typing import List

import pandas as pd

def createOrCleanDir(DIRPATH:str) -> None:
    """
    Create or clean a directory

    Args:
        DIRPATH (str): path to the direcotry which needs to be created or cleaned
    """
    if not os.path.exists(DIRPATH):
        print(f"{DIRPATH} does not exist. Creating ...")
        os.makedirs(DIRPATH)
    else:
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


def queryCSVFiles(DIRPATH:str) -> List[str]:
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
            if item.endswith(".csv"):
                fullPath = os.path.join(DIRPATH, item)
                result.append(fullPath)
    else:
        print(f"{DIRPATH} does not exist but is queried for CSV files")
    return result

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