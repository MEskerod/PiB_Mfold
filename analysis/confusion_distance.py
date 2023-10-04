import os, argparse
import matplotlib.pyplot as plt

def get_path_list(dir_name): 
    """
    Makes a list of all the files in the given directory. 
    Returns a list of "dir_name/file_name"
    """
    file_names = [os.path.join(dir_name, name) for name in os.listdir(dir_name)]
    return file_names

def read_dbn_file(file):
    """
    Read the dot bracket structure and sequence from a .dbn file. 
    The .dbn file has header lines that starts with #
    The first non-header line is the sequence and the second is the dot bracket structure

    Returns dotbracket structure, sequence
    """
    lines = []
    with open(file, 'r') as f: 
        for line in f: 
            if line.startswith('#'): 
                continue
            lines.append(line.strip())

    return lines[1], lines[0] #dot bracket, sequence

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

    return (len(structure1), F)

def distances(dir1, dir2): 
    
    structures1 = get_path_list(dir2)
    
    
    structures2 = get_path_list(dir1)
    
    
    lengths = []
    Fs = []

    for n in range(len(structures1)): 
        l, F = calculate_F(structures1[n], structures2[n])
        lengths.append(l)
        Fs.append(F)

    return lengths, Fs

