import os, time, sys

from general import make_dir, generate_random_sequence, calculate_slice_lengths, write_fasta

def time_consumption(func_name: str, func: callable, file_list: list[str]) -> list[float]:
    """
    Takes a function as input as uses that to fold all the  sequences from the files in file_list. 
    Returns a list with time consumption for folding each of the sequences

    Args: 
        - func_name: name of the function
        - func: the function to be used for folding structures
        - file_list: a list of fasta files that contains RNA sequences
    """
    print(f"\nFOLDING WITH {func_name.upper()}", file=sys.stderr)
    
    output_dir = make_dir(os.path.join("..", "structures", func_name + "_structures")) 
    
    times = []
    for index, file in enumerate(file_list): 
        print(f"Fold {file} - {index + 1}/{len(file_list)}", file=sys.stderr)
        start_time = time.time()
        func(file, output_dir)
        times.append(time.time() - start_time)
    return times

def run_functions_time(func_dict: dict, file_list: list[str]) -> dict: 
    """
    Takes a dictionary that contains functions for folding RNA, with the function name as key.
    For each function it folds the sequences in file_list
    Returns a dictionary with times for each function 

    Args: 
        - func_dict: dictionary of functions. Keys are function names and values are function 
        - file_list: a list of fasta files that contains RNA sequences
    """
    keys  = list(func_dict.keys())
    time_dict = {key:time_consumption(key, func_dict[key], file_list) for key in keys}
    return time_dict

def time_synthetic(n_sequences: int, min_length: int, max_length: int, func: callable) -> tuple[list[float], list[int]]:
    """
    To be able to compare the time for the time from the entire program sequence has to be read from a file and the db has to be written to a file
    Returns times for folding synthethic sequences and lengths of the sequences
    
    Args: 
        - n_sequences: The number of synthetic seuences that shoul be generated
        - min_length: Length of the shortest sequence that is generated 
        - max_length: Length of the longest seuence that is generated 
        - func: The function that is used to fold the synthetic sequence
    """
    
    print(f"\nFOLDING WITH SYNTHETIC SEQUENCES", file=sys.stderr)

    #Generate synthethic sequence and slices
    sequence = generate_random_sequence(max_length, ['A', 'C', 'G', 'U'])
    slices = calculate_slice_lengths(n_sequences, min_length, max_length)

    temporary_fasta = "temp.fasta"
    
    times = []
    lengths = []
    
    #Run algorithm and save time
    for index, n in enumerate(slices): 
        lengths.append(len(sequence[:n]))
        write_fasta(sequence[:n], "sythetic sequence", temporary_fasta)
        print(f"Fold sequence {index + 1}/{len(slices)}", file=sys.stderr)
        start_time = time.time()
        func("temp.fasta", "") 
        times.append(time.time() - start_time)
    
    #Remove tempoary files
    os.remove(temporary_fasta)
    os.remove("temp.dbn")

    return times, lengths
