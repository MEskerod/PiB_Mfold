import os, subprocess, time, random
import matplotlib.pyplot as plt
import pandas as pd

from confusion_distance import distances as calculate_f_distance


### FUNCTIONS FOR DIFFERENT PROGRAMS ###
"""
The functions have to take file and output directory as parameters. 
Each function runs the functionality of the program on 'file' and outputs the resulting dot bracket structure to a .txt file
"""

def run_Mfold_web(file, outdir):
    command = f"python3 cmd_Mfold_db.py {file} {outdir}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

def run_Mfold_orginal(file, outdir):
    outputfile = os.path.join(outdir, file)
    command = f"python3 Mfold/main.py -f {file} -lp 1988 -o {outputfile}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()


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

def write_txt(lines: list, outputfile): 
    """
    Writes lines in a txt file
    """
    with open(outputfile, 'w') as f: 
        for line in lines: 
            f.write(line + "\n")

def make_dir(dir_name):
    """
    """
    if not os.path.exists(dir_name):
            os.makedirs(dir_name)
    return dir_name

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

    temporary_txt = "temp.txt"
    output_dir = os.getcwd()
    
    times = []
    
    for n in slices: 
        write_txt([">synthetic sequence", sequence[:n]], temporary_txt) 
        start_time = time.time()
        #FIXME - Use the funtion we want to compare time between synthetic and real
        times.append(time.time() - start_time)
    os.remove(temporary_txt)
    return

### PLOT FUNCTIONS ###
def plot(lenghts: list, distances: dict): 
    #FIXME - Does not work and is just copied from different file
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

def distance_plot(): 
    return

def time_plot(): 
    return

### COLLECTING EVERYTHING INTO THE ACTUAL ANALYSIS ###
def main():
    ### SETTING UP ###
    file_list = get_path_list() #FIXME - Add directory

    ### GENERATE STRUCTURES WHERE NO TIME IS NEEDED ###

    ### TIME FOR MFOLD ###
    func_time_real = {}
    
    time_dict = run_functions_time(func_time_real, file_list) #Time for real data
    
    #TODO - With synthetic/random data 
    #NOTE - Which version of Mfold?
    #NOTE - Only Mfold or also compare time to other algorithms

    ### CALCULATE DISTANCES BETWEEN STRUCTURES ###
    #TODO - Compare my Mfold to true structures 
    #TODO - Compare other alorithms to true structures
    #NOTE - Do we want to compare my Mfold to other algorithms or just to true structures? 

    ### MAKE PLOTS ###
    #TODO - Make barplot of distances between structures
    #TODO - Make plot of time consumption of Mfold comparing true and synthetic structures
    #NOTE - Plot of comparing time between algorithms?


#if __name__ == '__main__': 
#    main()

time_synthetic(30, 50, 700)