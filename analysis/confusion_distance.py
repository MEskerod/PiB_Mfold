import os, argparse
import matplotlib.pyplot as plt
from itertools import combinations

from general import get_path_list, read_dbn_file

def distance(structure1, structure2): 
    true_positive = 0
    true_negative = 0
    false_positive = 0
    false_negative = 0
    
    if len(structure1) != len(structure2): 
        raise ValueError("Structures must have same length")

    for n in range(len(structure1)): 
        s1 = 1 if structure1[n] in '()' else 0
        s2 = 1 if structure2[n] in '()' else 0

        if s1 == 1 and s2 == 1: 
            true_positive += 1
        elif s1 == 1 and s2 == 0: 
            false_negative += 1
        elif s1 == 0 and s2 == 0: 
            true_negative += 1
        elif s1 == 0 and s2 == 1: 
            false_positive += 1
    
    return {"TP": true_positive, "FP": false_positive, "TN": true_negative, "FN": false_negative}

def F_score(CM: dict):
    precision = CM["TP"]/(CM["TP"]+CM["FP"])
    recall = CM["TP"]/(CM["TP"]+CM["FN"])

    F = 2 * (precision*recall)/(precision+recall)

    return F

def calculate_F(file1, file2): 
    structure1, _ = read_dbn_file(file1)
    structure2, _ = read_dbn_file(file2)

    confusion_matrix = distance(structure1, structure2)

    F = F_score(confusion_matrix)

    return F

def distances(dir1, dir2): 
    
    structures1 = get_path_list(dir2)
    
    
    structures2 = get_path_list(dir1)
    
    Fs = []

    for n in range(len(structures1)): 
        assert os.path.basename(structures1[n]) == os.path.basename(structures2[n])
        F = calculate_F(structures1[n], structures2[n])
        Fs.append(F)

    return Fs


def Fdistances(folders): 
    """
    """
    distance_dict = {}
    
    folder_combinations = list(combinations(folders, 2))

    for combination in folder_combinations: 
        folder1 = os.path.basename(combination[0]).split('_')[0]
        folder2 = os.path.basename(combination[1]).split('_')[0]
        name = folder1 + '_' + folder2
        distance = distances(combination[0], combination[1])
        distance_dict[name] = distance

    return distance_dict

