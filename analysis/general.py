import pandas as pd
import os, random

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord


### GENERAL FUNCTIONS ###
def get_path_list(dir_name): 
    """
    Makes a list of all the files in the given directory. 
    Returns a list of "dir_name/file_name"
    """
    file_names = [os.path.join(dir_name, name) for name in os.listdir(dir_name)]
    return file_names


def write_csv(data, output_file): 
    """
    Data is written to a .csv file
    """
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)

def write_fasta(sequence, seq_id, file_name):
    """
    """
    seq = Seq(sequence)
    record =  SeqRecord(seq, id = seq_id)

    with open(file_name, "w") as f: 
        SeqIO.write(record, f, "fasta")

def make_dir(dir_name):
    """
    """
    if not os.path.exists(dir_name):
            os.makedirs(dir_name)
    return dir_name



### ANALYSIS FUNCTIONS ###

def get_len_and_type(file_list): 
    """
    The files has to be named as ... (len_type_rest)
    """
    
    len_and_type = []

    for file in file_list: 
        splits = os.path.basename(file).split('_')
        len_and_type.append((splits[0], splits[1]))

    
    return len_and_type


def generate_random_sequence(length: int, alphabet: list): 
    """
    """
    random_sequence = "".join(random.choice(alphabet) for _ in range(length))
    return random_sequence

def create_quadratic_function(y_intercept: int, point: tuple): 
    """
    Returns the quadratic function. 
    The quadratic function that is returned will have a positive a value and will be centered around 0 (b=0)

    f(x) = a*x^2 + b*x + c and since b = 0, it is just f(x) = a*x^2 + b
    """
    c = y_intercept

    x1, y1 = point[0], point[1]

    a = (y1 - c)/x1**2

    def quadratic(x): 
        return a * x**2 + c
    
    return quadratic

def calculate_slice_lengths(num_slices, min_length, initial_length):
    slice_lengths = []
    quadratic = create_quadratic_function(min_length, (num_slices, initial_length))
    for x in range(1, num_slices+1): 
        slice_lengths.append(int(quadratic(x)))
    return slice_lengths