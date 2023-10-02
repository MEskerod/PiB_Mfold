import time, os, argparse
import matplotlib.pyplot as plt
import pandas as pd

from Mfold import(read_parameters, 
                  read_fasta, 
                  prepare_input,
                  fold_rna,
                  backtrack)

def get_path_list(dir_name): 
    """
    Makes a list of all the files in the given directory. 
    Returns a list of "dir_name/file_name"
    """
    file_names = [os.path.join(dir_name, name) for name in os.listdir(dir_name)]
    return file_names

def read_file(file_path):
    """
    Opens a fasta file and reads and prepares the sequence. 
    An error is given if the sequence contains not allowed letters 
    """
    with open(file_path, 'r') as f: 
        sequence, name = read_fasta(f)
    sequence = prepare_input(sequence)
    return sequence, name

def complete_fold(path): 
    """
    Add all the functionality to fold a RNA into RNA secondary structure using Mfold
    Returns the length of the sequence and the fold
    """
    sequence, name =  read_file(path)
    N = len(sequence)
    parameters = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    W, V = fold_rna(sequence, parameters)
    fold = backtrack(W, V, parameters, sequence)
    return N, fold, name

def write_fold_to_file(output_dir, file, fold, name): 
    """
    Takes the dot bracket structure and writes it to a .txt file
    The .txt file has the same name as the input file and is placed in a folder called 'structures' in output_dir
    """
    output_file = os.path.join(output_dir, "structures", os.path.splitext(os.path.basename(file))[0]+".txt")
    with open(output_file, 'w') as f: 
        f.write(f">{name}\n")
        f.write(f"{fold}\n")

def running_times(files, output_path): 
    """
    For each file in the input directory the running time for folding the sequence is found. 
    The resulting fold is written to a .txt file
    A dictionary is returned containing the lengths and corresponding running times of the sequences. 
    """
    times = []
    lengths = []

    for file in files: 
        print(f"Folding {file}")
        starttime = time.time()
        N, fold, name = complete_fold(file)
        endtime = time.time()
        times.append(endtime-starttime)
        lengths.append(N)
        write_fold_to_file(output_path, file, fold, name)

    data = {'length': lengths, 
            'time': times}
    
    return data

def plot_times(data: dict, output_path):
    """
    Running time is plotted as a function of sequence length
    """
    x = data["length"] 
    y = data["time"]

    plt.scatter(x, y, color = 'black')
    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.grid(True)
    plt.savefig(os.path.join(output_path, "runningtimes.jpeg"))

def write_csv(data, output_file): 
    """
    Data is written to a .csv file
    """
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)


def main(): 
    """
    Calculates and plots the running time for all sequences in files in input_dir. 
    In output_dir are the folds of all the sequences, a .csv file with the times and a plot of the running time
    """
     #Setting up the option parsing using the argparse module
    argparser = argparse.ArgumentParser(
        description="" )
    #Adding arguments
    #TODO - Add description/help for command line options
    argparser.add_argument('input_dir') 
    argparser.add_argument('output_dir')

    args = argparser.parse_args()
    
    output_dir = args.output_dir
    
    files = get_path_list(args.input_dir)
    data = running_times(files, output_dir)
    plot_times(data, output_dir)
    write_csv(data, os.path.join(output_dir, "times.csv"))





if __name__ == '__main__': 
    main() 