import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
import os
import numpy as np

### PLOT FUNCTIONS ###
def plot_synthetic_real_times(synthethic_lengths, real_lengths, synthetic_time, real_time, output_path):
    """
    Running time is plotted as a function of sequence length, for synthethic sequences vs real sequences
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

def plot_Nussinov(lengths, time, output_path): 
    """
    Running time is plotted as a function of sequence length for Nussinov algorithm
    """
    plt.figure()
    plt.scatter(lengths, time, s=20, edgecolors='blue', facecolors = 'none', linewidths=1)
    plt.plot(lengths, time, color = 'blue', linestyle = '--', linewidth = 0.8)
    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.grid(True, linestyle = '--')
    plt.savefig(os.path.join(output_path, "Nussinov_times.jpeg"))
    return

def plot_MfoldOriginal_newest(lengths, original_time, newest_time, output_path): 
    """
    Running time is plotted as a function of sequence length for newest implemented version of Mfold compared to the original
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

def plot_distances(distance_dict: dict, key_label, output_file):
    """
    Distances are plotted as bar plot. 
    Distances for same sequence are groupped together. 
    The sequences are sorted according to type and then length
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

def plot_single_distance(): 
    return