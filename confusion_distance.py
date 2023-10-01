import os, argparse
import matplotlib.pyplot as plt

def get_path_list(dir_name): 
    """
    Makes a list of all the files in the given directory. 
    Returns a list of "dir_name/file_name"
    """
    file_names = [os.path.join(dir_name, name) for name in os.listdir(dir_name)]
    return file_names

def read_dot_bracket(file):
    """
    Read the dot bracket structure in the text file into a string
    """
    with open(file, 'r') as f: 
        file_content = f.readlines()
    
    dotbracket = "".join([line.strip() for line in file_content if not line.startswith(">")])
    return dotbracket

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
    structure1 = read_dot_bracket(file1)
    structure2 = read_dot_bracket(file2)

    confusion_matrix = distance(structure1, structure2)

    F = F_score(confusion_matrix)

    return (len(structure1), F)

def plot(lenghts: list, distances: dict): 
    colors = ["red", "blue", "green", "yellow"]
    keys = list(distances.keys())

    x = [str(l) for l in lenghts]

    for n in range(len(keys)):
        key = keys[n]
        plt.bar(x, distances[key], color=colors[n], label = key)
    plt.ylabel("f score")
    plt.xlabel("Sequence length")
    plt.grid(True)
    plt.legend()
    plt.savefig("distances.jpeg")
    return

def main(true_dir, my_dir): 
    
    my_structures = get_path_list(my_dir)
    
    
    true_structures = get_path_list(true_dir)
    
    
    lengths = []
    Fs = []

    for n in range(len(my_structures)): 
        l, F = calculate_F(my_structures[n], true_structures[n])
        lengths.append(l)
        Fs.append(F)
    
    distances = {"My/True": Fs}

    plot(lengths, distances)

    return lengths, Fs


path1 = "true_structures"
path2 = "results/structures"

print(main(path1, path2))
