import matplotlib.pyplot as plt

from general import get_path_list

from time_functions import run_functions_time

from confusion_distance import distances as calculate_f_distance


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