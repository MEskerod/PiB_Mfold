from general import get_path_list, write_csv, get_len_and_type, read_csv
from running_fold import run_Mfold_web, run_Mfold_orginal, run_Mfold_newest, run_Nussinov
from time_functions import run_functions_time, time_synthetic
from confusion_distance import Fdistances
from check_files import check_files
from plots import plot_synthetic_real_times, plot_Nussinov, plot_MfoldOriginal_newest, plot_distances

import os, sys


### COLLECTING EVERYTHING INTO THE ACTUAL ANALYSIS ###
def main() -> None:
    """
    Performs (almost) all the analysis for Mfold PiB project
    """
    ### SETTING UP ###
    file_list = get_path_list("../sequences")
   
    os.makedirs("../results", exist_ok=True)

    #Get lenght and type for seuences
    lengths, types = get_len_and_type(file_list)

    ### TIME FOR ALGORITHMS ###
    #Generate time for all and save to .csv. Plotd an be generated or changed later
    func_time_real = {"OriginalMfold": run_Mfold_orginal, 
                      "NewestMfold": run_Mfold_newest,  
                      "Nussinov": run_Nussinov, 
                      "MfoldWebversion": run_Mfold_web}

    
    #Generate times for all functions
    time_dict = run_functions_time(func_time_real, file_list) #Time for real data
    time_dict["Length"] = lengths

    #Save times to .csv
    write_csv(time_dict, "../results/time_table.csv")

    #Find time for synthetic seuences and add to time_dict 
    synthetic_time, synthetic_lengths = time_synthetic(40, 60, 800, run_Mfold_newest)
    synthetic_dict = {"NewestMfoldSynthetic": synthetic_time, "Length": synthetic_lengths}
    write_csv(synthetic_dict, "../results/syntethic_times.csv")

    """
    #FOR RE-LOADING DICTIONARIES
    synthetic_dict = read_csv("../results/syntethic_times.csv")
    time_dict = read_csv("../results/time_table.csv")
    """

    ### CALCULATE DISTANCES BETWEEN STRUCTURES ###
    #Get a list of all the folders with structures: 
    folders = [os.path.join("../structures", subdir) for subdir in os.listdir("../structures")]

    #Check content 
    print("\nCHECKING FILES", file=sys.stderr)
    for folder in folders: 
        check_files("../sequences", folder) 


    #Compare all to each other and save to .csv to be used later if needed
    print("\nCALCULATE DISTANCES", file=sys.stderr)
    distance_dict = Fdistances(folders)
    distance_dict["length"] = lengths
    distance_dict ["type"] = types
    write_csv(distance_dict, "../results/distance_table.csv")

    ### MAKE PLOTS ###
    print("\nMAKE PLOTS", file=sys.stderr)
    plot_synthetic_real_times(synthetic_dict["Length"], time_dict["Length"], synthetic_dict["NewestMfoldSynthetic"], time_dict["NewestMfold"], "../results") 
    plot_MfoldOriginal_newest(time_dict["Length"], time_dict["OriginalMfold"], time_dict["NewestMfold"], "../results")
    plot_Nussinov(time_dict["Length"], time_dict["Nussinov"], "../results")

    distances_to_plot1 = [("true_NewestMfold", "Newest version"), ("true_OriginalMfold", "Original version"), ("true_Nussinov", "Nussinov"), ("true_MfoldWebversion", "Mfold web")]
    plot_distance_dict1 = {key[0]: distance_dict[key[0]] for key in distances_to_plot1}
    plot_distance_dict1["length"] = lengths
    plot_distance_dict1["type"] = types

    distances_to_plot2 = [("NewestMfold_OriginalMfold", "Original version"), ("NewestMfold_Nussinov", "Nussinov"), ("NewestMfold_MfoldWebversion", "Mfold web")]
    plot_distance_dict2 = {key[0]: distance_dict[key[0]] for key in distances_to_plot2}
    plot_distance_dict2["length"] = lengths
    plot_distance_dict2["type"] = types

    plot_distances(plot_distance_dict1, distances_to_plot1, "../results/distances_true.jpeg")
    plot_distances(plot_distance_dict2, distances_to_plot2, "../results/distances_newest.jpeg")


if __name__ == '__main__': 
    main()
