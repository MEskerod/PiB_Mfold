import os, time

from general import make_dir, generate_random_sequence, calculate_slice_lengths, write_fasta

def time_consumption(func_name, func, file_list):
    """
    """
    output_dir = make_dir(func_name + "_structures")
    
    times = []
    for file in file_list: 
        start_time = time.time()
        func(file, output_dir)
        times.append(time.time() - start_time)
    return times

def run_functions_time(func_dict: dict, file_list: list): 
    """
    """
    keys  = list(func_dict.keys())
    time_dict = {key:time_consumption(key, func_dict[key], file_list) for key in keys}
    return time_dict

def time_synthetic(n_sequences, min_length, max_length):
    """
    To be able to compare the time for the time from the entire program sequence has to be read from a file and the db has to be written to a file
    """
    
    sequence = generate_random_sequence(max_length, ['A', 'C', 'G', 'U'])
    slices = calculate_slice_lengths(n_sequences, min_length, max_length)

    temporary_txt = "temp.fasta"
    output_dir = os.getcwd()
    
    times = []
    
    for n in slices: 
        write_fasta(sequence[:n], "sythetic sequence", temporary_txt)
        start_time = time.time()
        #FIXME - Use the funtion we want to compare time between synthetic and real
        times.append(time.time() - start_time)
    os.remove(temporary_txt)
    return
