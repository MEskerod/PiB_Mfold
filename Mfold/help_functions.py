import ast, argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import TextIOWrapper

#TODO - Read and update all descriptons to match the content!
#TODO - Add types for all parameters

### HELPER FUNCTIONS USED IN MAIN ###
def read_fasta(input) -> str:
    """
    Reads in a FASTA-file and returns the sequence
    If there is more than one sequence in the FASTA file it gives an error
    """
    records = list(SeqIO.parse(input, 'fasta'))
    
    if len(records) > 1: 
        raise ValueError("FASTA file contains more than one sequence")

    return str(records[0].seq), records[0].id

def prepare_input(input: str) -> str: 
    """
    Strips new line in the end
    Makes sure that all letters are in upper case
    A value error is raised if any invalid letters are found in the sequence.  
    """
    allowed = set(['A', 'C', 'G', 'U', 'N'])
    input = input.strip().upper().replace('T', 'U')

    if any(letter not in allowed for letter in input): 
        raise ValueError("Invalid letters found in the sequence")
    
    return input

def read_general_parameters(file_loop, file_stacking): 
    """
    Read parameters from to .csv files:
    - One containing the loop parameters 
    - A second one containg the parameters for different types og base pairing 
    The .csv files are converted to Pandas tables and returned
    """
    try:
        loops = pd.read_csv(file_loop) 
        stacking = pd.read_csv(file_stacking, index_col=0)
    except FileNotFoundError:
        raise FileNotFoundError("One or both parameter files not found")
    except pd.errors.EmptyDataError:
        raise ValueError("One or both parameter files are empty or in an unexpected file format")
    return loops, stacking

def write_dbn(name, sequence, fold, outfile: TextIOWrapper):
    """
    """
    outfile.write(f"#Name: {name}\n")
    outfile.write(f"#Length: {len(sequence)}\n")
    outfile.write(sequence + "\n")
    outfile.write(fold + "\n")


def parse_asymmetry_parameters(input_string: str): 
    """
    """
    try: 
        #Split into the f and penalty max parts
        arguments = input_string.split(';') 

        #Evaluate the input safely and evaluate literals. Can hande lists. 
        f = ast.literal_eval(arguments[0])

        penalty_max = int(arguments[1])

        return f, penalty_max
    
    except Exception as e: 
        raise argparse.ArgumentTypeError(f"Invalid input: {input_string}. Error: {str(e)}")

