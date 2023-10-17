from general import get_path_list, write_csv, get_len_and_type
from running_fold import run_Mfold_web, run_Mfold_orginal, run_Mfold_newest, run_Nussinov
from time_functions import run_functions_time, time_synthetic
from confusion_distance import Fdistances
from check_files import check_files
from plots import plot_synthetic_real_times, plot_Nussinov, plot_MfoldOriginal_newest, plot_distances

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
    #TODO - Change distance to take care of other types of brackets!
    distance_dict = Fdistances(folders)
    distance_dict["length"] = lengths
    distance_dict ["type"] = types
    write_csv(distance_dict, "../results/distance_table.csv")

    ### MAKE PLOTS ###
    #plot_synthetic_real_times(lengths, time_dict["NewestMfoldSynthetic"], time_dict["NewestMfold"], "../results") #NOTE - Change times!!!
    plot_MfoldOriginal_newest(lengths, time_dict["OriginalMfold"], time_dict["NewestMfold"], "../results")
    plot_Nussinov(lengths, time_dict["Nussinov"], "../results")

    distances_to_plot1 = [("NewestMfold_true", "Newest version"), ("OriginalMfold_true", "Original version"), ("MfoldWebversion_true", "Mfold web"), ("Nussinov_true", "Nussinov")]
    plot_distance_dict1 = {key[0]: distance_dict[key[0]] for key in distances_to_plot1}
    plot_distance_dict1["length"] = lengths
    plot_distance_dict1["type"] = types

    distances_to_plot2 = [("NewestMfold_OriginalMfold", "Original version"), ("MfoldWebversion_NewestMfold", "Mfold web"), ("NewestMfold_Nussinov", "Nussinov")]
    plot_distance_dict2 = {key[0]: distance_dict[key[0]] for key in distances_to_plot2}
    plot_distance_dict2["length"] = lengths
    plot_distance_dict2["type"] = types

    plot_distances(plot_distance_dict1, distances_to_plot1, "../results/distances_true.jpeg")
    plot_distances(plot_distance_dict2, distances_to_plot2, "../results/distances_newest.jpeg")
    
    #TODO - Find out if we want more distance comparisons

    #FIXME - Before running!!! Change paths! And check for paths in other files! (time_functions)


if __name__ == '__main__': 
    main()
