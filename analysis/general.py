import pandas as pd
import os, random

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord


### GENERAL FUNCTIONS ###
def get_path_list(dir_name: str) -> list[str]: 
    """
    Makes a list of all the files in the given directory. 
    Returns a list of "dir_name/file_name"
    """
    file_names = [os.path.join(dir_name, name) for name in sorted(os.listdir(dir_name), key=lambda x:(int(x.split('_')[0]), x.split('_')[1], x.split('_')[2]))]
    return file_names


def write_csv(data, output_file: str) -> None: 
    """
    Data is written to a .csv file
    The input data can be a dictionary, with the keys being the column names and values being the columns
    """
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)

def read_csv(input_file: str) -> dict: 
    """
    Reads a .csv file and returns the data as a dictionary
    """
    df = pd.read_csv(input_file)
    data_dict = df.to_dict(orient='list')
    return data_dict

def write_fasta(sequence: str, seq_id: str, file_name: str) -> None:
    """
    Writes a sequence to a fasta file

    Args:
        - sequence: the sequence that has to be written to a file 
        - seq_id: the id/name of the sequence
        - file_name: the name of the file. Has to have the extension .fasta
    """
    seq = Seq(sequence)
    record =  SeqRecord(seq, id = seq_id)

    with open(file_name, "w") as f: 
        SeqIO.write(record, f, "fasta")

def make_dir(dir_name: str) -> str:
    """
    Make a directory, if the directory does not exist
    """
    if not os.path.exists(dir_name):
            os.makedirs(dir_name)
    return dir_name

def read_dbn_file(file: str) -> tuple[str, str]:
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

def read_fasta(input: str) -> str:
    """
    Reads in a FASTA-file and returns the sequence
    If there is more than one sequence in the FASTA file it gives an error
    """
    records = list(SeqIO.parse(input, 'fasta'))
    
    if len(records) > 1: 
        raise ValueError("FASTA file contains more than one sequence")

    return str(records[0].seq), records[0].id

### ANALYSIS FUNCTIONS ###

def get_len_and_type(file_list: list[str]) -> tuple[list[int], list[str]]: 
    """
    Takes a list of files and returns two list. One containing the lengths of the sequences in the file an a second containing the type of RNA. 
    The files has to be names as some/path/lenght_type_rest
    """
    
    lenghts = []
    types = []

    for file in file_list: 
        splits = os.path.basename(file).split('_')
        lenghts.append(int(splits[0]))
        types.append(splits[1])
        
    return lenghts, types


def generate_random_sequence(length: int, alphabet: list) -> str: 
    """
    Generates a random seuence of a given length. 
    
    Args:
        - length: length of the sequence that is generated
        - alphabet: the characters from which to generate the sequence
    """
    random_sequence = "".join(random.choice(alphabet) for _ in range(length))
    return random_sequence

def create_quadratic_function(y_intercept: int, point: tuple[int, int]) -> callable: 
    """
    Returns the quadratic function. 
    The quadratic function that is returned will have a positive a value and will be centered around 0 (b=0)

    f(x) = a*x^2 + b*x + c and since b = 0, it is just f(x) = a*x^2 + b

    Args: 
        - y_intercept: the point where the quadratic funtion intercepts the y-axis
        - point: a point that fullfils = f(x) = y
    """
    c = y_intercept

    x1, y1 = point[0], point[1]

    a = (y1 - c)/x1**2

    def quadratic(x): 
        return a * x**2 + c
    
    return quadratic

def calculate_slice_lengths(num_slices: int, min_length: int, initial_length: int) -> list[int]:
    """
    Calculate the lengths of the slices, to obtain a given number of slices of lengths between a minimum and maximum length, spaced according to a quadratic function. 

    Args:
        - num_slices: Number of slices wanted 
        - min_length: Length of the shortest slice 
        - initial_length: Length of the sequence to slice from, which is also equal to the maximum length
    """
    slice_lengths = []
    quadratic = create_quadratic_function(min_length, (num_slices, initial_length))
    for x in range(1, num_slices+1): 
        slice_lengths.append(int(quadratic(x)))
    return slice_lengths