from general import get_path_list, write_csv
from running_fold import run_Mfold_web, run_Mfold_orginal, run_Mfold_newest, run_Nussinov
from time_functions import run_functions_time, time_synthetic
from confusion_distance import distances as calculate_f_distance

import os


### COLLECTING EVERYTHING INTO THE ACTUAL ANALYSIS ###
def main():
    ### SETTING UP ###
    file_list = get_path_list("../sequences") #NOTE - Change for test!!!
    os.makedirs("../results", exist_ok=True)

    ### TIME FOR ALGORITHMS ###
    #Generate time for all and save to .csv. Plotd an be generated or changed later
    func_time_real = {"Original Mfold": run_Mfold_orginal, 
                      "Newest Mfold": run_Mfold_newest, 
                      "Mfold web version": run_Mfold_web, 
                      "Nussinov": run_Nussinov}
    
    #Generate times for all functions
    time_dict = run_functions_time(func_time_real, file_list) #Time for real data

    #Find time for synthetic seuences and add to time_dict
    synthetic_time = time_synthetic(50, 50, 800, run_Mfold_newest)
    time_dict["Newest Mfold synthetic"] = synthetic_time

    #Save times to .csv
    write_csv(time_dict, "../results/time_table.csv")


    ### CALCULATE DISTANCES BETWEEN STRUCTURES ###
    #Compare all to each other and save to .csv to be used later if needed
    
    #TODO - Compare my Mfold to true structures 
    #TODO - Compare other algorithms to true structures
    #NOTE - Do we want to compare my Mfold to other algorithms or just to true structures? 

    ### MAKE PLOTS ###
    #TODO - Make barplot of distances between structures
    #TODO - Make plot of time consumption of Mfold comparing true and synthetic structures
    #NOTE - Plot of comparing time between algorithms?


if __name__ == '__main__': 
    main()
