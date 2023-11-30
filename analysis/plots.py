import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import os
import numpy as np

### PLOT FUNCTIONS ###
def plot_synthetic_real_times(synthethic_lengths: list[int], real_lengths: list[int], synthetic_time: list[float], real_time: list[float], output_path: str) -> None:
    """
    Running time is plotted as a function of sequence length, for synthethic sequences vs real sequences

    Args: 
        - synthethic_lengths: list containing lengths of synthethic sequences (same order as syntethic_time) 
        - real_lengths: list containing lengths of real seqeuences (same order as real_time)
        - synthethic_time: list containing time used for folding synthethic sequences (same order as synthethic_lengths) 
        - real_time: list containing time used for folding real sequences (same order as real_lengths)
        - output_path: path to where the plot is saved
    """
    plt.figure()
    handles = [] #For adding descriptions to legend (lines wth markers)
    
    plt.scatter(real_lengths, real_time, s=20, edgecolors='blue', facecolors = 'none', linewidths=1)
    plt.plot(real_lengths, real_time, color = 'blue', linestyle = '--', linewidth = 0.8)
    handles.append(Line2D([0], [0], color = 'blue', linestyle='--', linewidth=0.8, marker='o', markerfacecolor='none', markeredgecolor='blue', label = "Real RNA sequences"))

    plt.scatter(synthethic_lengths, synthetic_time, s=20, edgecolors= 'red', facecolors = 'none', linewidths=1)
    plt.plot(synthethic_lengths, synthetic_time, color = 'red', linestyle = '--', linewidth = 0.8)
    handles.append(Line2D([0], [0], color='red', linestyle='--', linewidth=0.8, marker='o', markerfacecolor='none', markeredgecolor='red', label = "Synthetic sequences"))

    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.legend(handles=handles)
    plt.grid(True, linestyle = '--')
    plt.savefig(os.path.join(output_path, "synthetic_vs_real_times.jpeg"))

def plot_Nussinov(lengths: list[int], time: list[float], output_path: str) -> None: 
    """
    Running time is plotted as a function of sequence length for Nussinov algorithm

    Args: 
        - lenghts: list containing lengths of the sequences (same order as time)
        - time: list containing times for folding the sequences (same order as lengths)
        - output_path: path to where the plot should be saved
    """
    plt.figure()
    plt.scatter(lengths, time, s=20, edgecolors='blue', facecolors = 'none', linewidths=1)
    plt.plot(lengths, time, color = 'blue', linestyle = '--', linewidth = 0.8)
    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.grid(True, linestyle = '--')
    plt.savefig(os.path.join(output_path, "Nussinov_times.jpeg"))
    return

def plot_MfoldOriginal_newest(lengths: list[int], original_time: list[float], newest_time: list[float], output_path: str): 
    """
    Running time is plotted as a function of sequence length for newest implemented version of Mfold compared to the original

    Args: 
        - lenghts: list containing lengths of the sequences
        - original_time: list containing times for folding the sequences with original verion of Mfold (same order as lengths)
        - newest_time: list containing times for folding the sequences with newest version of Mfold (same order as lengths)
        - output_path: path to where the plot should be saved
    """
    plt.figure()
    handles = [] #For adding descriptions to legend (lines wth markers)

    plt.scatter(lengths, original_time, s=20, edgecolors='blue', facecolors = 'none', linewidths=1)
    plt.plot(lengths, original_time, color = 'blue', linestyle = '--', linewidth = 0.8)
    handles.append(Line2D([0], [0], color = 'blue', linestyle='--', linewidth=0.8, marker='o', markerfacecolor='none', markeredgecolor='blue', label = "Original Mfold"))

    plt.scatter(lengths, newest_time, s=20, edgecolors= 'red', facecolors = 'none', linewidths=1)
    plt.plot(lengths, newest_time, color = 'red', linestyle = '--', linewidth = 0.8)
    handles.append(Line2D([0], [0], color='red', linestyle='--', linewidth=0.8, marker='o', markerfacecolor='none', markeredgecolor='red', label = "Mfold with improvements"))

    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.legend(handles = handles)
    plt.grid(True, linestyle = '--')
    plt.savefig(os.path.join(output_path, "original_vs_new_times.jpeg"))
    return

def plot_distances(distance_dict: dict, key_label: list[tuple], output_file: str) -> None:
    """
    Distances are plotted as bar plot. 
    Distances for same sequence are groupped together. 
    The sequences are sorted according to type and then length

    args: 
        - distance_dict: a dictionary containing the distances. key is name of the distance pair and values a lists containing the calculated distances
        - key_label: a list containing tuples, where the first element in the tuple is the name of the distance in distance_dict and second element is the one used for the plot
        - output_file: name (and path) of the output file
    """
    #Format data
    data = [{key: distance_dict[key][n] for key in list(distance_dict.keys())} for n in range(len(distance_dict["length"]))]

    def custom_sort_key(item):
        return (item["type"] == 'other', item["type"], item["length"]) #To push the type other to the end

    
    #Sort the data by type then by length
    data.sort(key=custom_sort_key)

    types = [d["type"] for d in data]
    lengths = [d["length"] for d in data]
    distances = {key[1]: [d[key[0]] for d in data] for key in key_label}

    tick_labels = [f"{types[n]}, {lengths[n]}" for n in range(len(data))]

    #Get unique RNA types and their position
    unique_types, type_positions = np.unique(types, return_index=True)
    
    x = np.arange(len(data)) #Make locations for labels
    width = 6.61/len(data)

    multiplier = 0

    fig, ax = plt.subplots(constrained_layout=True, figsize = (25, 12))

    #Plot groups
    for label, distance in distances.items(): 
        offset = width * multiplier
        ax.bar(x + offset, distance, width, label = label, edgecolor = 'black', align='edge')
        multiplier += 1

    #Add vertical lines
    for position in type_positions[1:]: 
        ax.axvline(x=position, linestyle = '--')

    #Set xticks and rotate
    ax.set_xticks(x + width * (multiplier-1)/2, tick_labels)
    ax.set_xticklabels(tick_labels, rotation = 90, fontsize = 16)

    ax.tick_params(axis='y', labelsize=18)
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}")) #Round y ticks to to decimal
    
    ax.legend(loc = 'upper left', bbox_to_anchor = (1.01, 1), fontsize = 16) #Set position of legend
    ax.set_ylabel("f1 score", fontsize = 22)
    
    #Add horizontal grid lines behind bars
    ax.set_axisbelow(True)
    ax.grid(axis='y', linestyle = '--') 

    plt.savefig(output_file)