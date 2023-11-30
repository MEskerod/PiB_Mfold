import os
from itertools import permutations

from general import get_path_list, read_dbn_file

def distance(structure1: str, structure2 : str) -> dict: 
    """
    Calculates the f1 score between structure1 and structure2, with structure1 as the ground truth
    Structure can contain the elements (, ), [, ], {, }, <, > and

    Returns a dictionary containing the number of false positives and negatives and true negatives and positives
    """
    allowed = {'(', ')', '[', ']', '{', '}', '<', '>', '.'}
    brackets = {'(', ')', '[', ']', '{', '}', '<', '>'}
    
    true_positive = 0
    true_negative = 0
    false_positive = 0
    false_negative = 0
    
    if len(structure1) != len(structure2): 
        raise ValueError("Structures must have same length")

    for n in range(len(structure1)): 
        if structure1[n] not in allowed or structure2[n] not in allowed: 
            raise ValueError(f"ERROR: Prohibited chrachters in one or both sequences at position {n}")
        
        #Find out if charachters is positives or negatives. 
        #Brackets/paired bases are positive and yields a 1
        #Dots/unpaired bases are negative and yields a 0
        s1 = 1 if structure1[n] in brackets else 0
        s2 = 1 if structure2[n] in brackets else 0

        if s1 == 1 and s2 == 1: 
            true_positive += 1
        elif s1 == 1 and s2 == 0: 
            false_negative += 1
        elif s1 == 0 and s2 == 0: 
            true_negative += 1
        elif s1 == 0 and s2 == 1: 
            false_positive += 1
    
    return {"TP": true_positive, "FP": false_positive, "TN": true_negative, "FN": false_negative}

def f_score(CM: dict, epsilon: float=1e-7) -> float:
    """
    Takes a dictionary containing the count of false positives and negatives and true negatives and positives
    epsilon is added to avoid division by zero

    args: 
        - CM: dictionary with TP = true positive, FP = false positive, FP = false positive, FN = false negative
        - epsilon: a small number that is added to avoid division by zero (default = 1e-7)
    """
    precision = CM["TP"]/(CM["TP"]+CM["FP"]+epsilon)
    recall = CM["TP"]/(CM["TP"]+CM["FN"]+epsilon)

    f1 = 2 * (precision*recall)/(precision+recall+epsilon)

    return f1

def calculate_f1(file1: str, file2: str) -> float: 
    """
    Calculates the f1 score between the structures in file1 and file2. 
    Files must be dbn files with header lines denoted by #
    """
    structure1, _ = read_dbn_file(file1)
    structure2, _ = read_dbn_file(file2)

    confusion_matrix = distance(structure1, structure2)

    f = f_score(confusion_matrix)

    return f

def distances(dir1: str, dir2: str) -> list[float]: 
    """
    Takes to directories and calculates the f1 score between the files in them. 
    The direcotries must contain the same structures and with names that ensures that get_path_list returns them in the same order 
    (Prefeably named the same)
    """
    #Get lists of files with structures
    structures1 = get_path_list(dir1)
    structures2 = get_path_list(dir2)
    
    Fs = []

    for n in range(len(structures1)): 
        assert os.path.basename(structures1[n]) == os.path.basename(structures2[n])
        F = calculate_f1(structures1[n], structures2[n])
        Fs.append(F)

    return Fs


def Fdistances(folders: list[str]) -> dict: 
    """
    Takes a list of directories. 
    It returns a dictionary containing the f1 scores between the sequences in all posible combinations of the directories in 'folders'
    """
    distance_dict = {}
    
    folder_combinations = list(permutations(folders, 2)) #Permutations will give all possible combination

    for combination in folder_combinations: 
        folder1 = os.path.basename(combination[0]).split('_')[0]
        folder2 = os.path.basename(combination[1]).split('_')[0]
        name = folder1 + '_' + folder2
        distance = distances(combination[0], combination[1])
        distance_dict[name] = distance

    return distance_dict

