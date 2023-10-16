from general import get_path_list, write_csv, get_len_and_type
from running_fold import run_Mfold_web, run_Mfold_orginal, run_Mfold_newest, run_Nussinov
from time_functions import run_functions_time, time_synthetic
from confusion_distance import Fdistances
from check_files import check_files
from plots import plot_synthetic_real_times

import os


### COLLECTING EVERYTHING INTO THE ACTUAL ANALYSIS ###
def main():
    ### SETTING UP ###
    file_list = get_path_list("../examples") #../sequences #NOTE - Change
    os.makedirs("../results", exist_ok=True)

    #Get lenght and type for seuences
    lengths, types = get_len_and_type(file_list)

    ### TIME FOR ALGORITHMS ###
    #Generate time for all and save to .csv. Plotd an be generated or changed later
    func_time_real = {"OriginalMfold": run_Mfold_orginal, 
                      "NewestMfold": run_Mfold_newest, 
                      "MfoldWebversion": run_Mfold_web, 
                      "Nussinov": run_Nussinov}
    
    
    #Generate times for all functions
    time_dict = run_functions_time(func_time_real, file_list) #Time for real data

    #Find time for synthetic seuences and add to time_dict #NOTE - Uncomment
    #synthetic_time = time_synthetic(50, 50, 800, run_Mfold_newest)
    #time_dict["NewestMfoldSynthetic"] = synthetic_time

    time_dict["Length"] = lengths

    #Save times to .csv
    write_csv(time_dict, "../results/time_table.csv")

    ### CALCULATE DISTANCES BETWEEN STRUCTURES ###
    #Get a list of all the folders with structures: 
    folders = [os.path.join("../structures_example", subdir) for subdir in os.listdir("../structures_example")] #../structures #NOTE - Change

    #Check content 
    for folder in folders: 
        check_files("../examples", folder) #../sequences #NOTE - Change


    #Compare all to each other and save to .csv to be used later if needed
    distance_dict = Fdistances(folders)
    distance_dict["length"] = lengths
    distance_dict ["type"] = types
    write_csv(distance_dict, "../results/distance_table.csv")

    ### MAKE PLOTS ###
    plot_synthetic_real_times(lengths, time_dict["NewestMfold"], time_dict["OriginalMfold"], "../results") #NOTE - Change times!!!


    #TODO - Make barplot of distances between structures
    #NOTE - Plot of comparing time between algorithms?


if __name__ == '__main__': 
    main()
