import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import numpy as np

### PLOT FUNCTIONS ###
def plot_synthetic_real_times(lengths, synthetic_time, real_time, output_path):
    """
    Running time is plotted as a function of sequence length
    """
    plt.figure()
    plt.scatter(lengths, real_time, label = "Real RNA sequences", s=20, edgecolors='blue', facecolors = 'none', linewidths=1)
    plt.plot(lengths, real_time, color = 'blue', linestyle = '--', linewidth = 0.8)
    plt.scatter(lengths, synthetic_time, label = "Synthetic sequences", s=20, edgecolors= 'red', facecolors = 'none', linewidths=1)
    plt.plot(lengths, synthetic_time, color = 'red', linestyle = '--', linewidth = 0.8)
    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.legend()
    plt.grid(True, linestyle = '--')
    plt.savefig(os.path.join(output_path, "synthetic_vs_real_times.jpeg"))

def plot_Nussinov(lengths, time, output_path): 
    """
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
    """
    plt.figure()
    plt.scatter(lengths, original_time, label = "Original Mfold", s=20, edgecolors='blue', facecolors = 'none', linewidths=1)
    plt.plot(lengths, original_time, color = 'blue', linestyle = '--', linewidth = 0.8)
    plt.scatter(lengths, newest_time, label = "Mfold with improvements", s=20, edgecolors= 'red', facecolors = 'none', linewidths=1)
    plt.plot(lengths, newest_time, color = 'red', linestyle = '--', linewidth = 0.8)
    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.legend()
    plt.grid(True, linestyle = '--')
    plt.savefig(os.path.join(output_path, "original_vs_new_times.jpeg"))
    return

def plot_distances(distance_dict: dict, key_label, output_file):
    #Format data
    data = [{key: distance_dict[key][n] for key in list(distance_dict.keys())} for n in range(len(distance_dict["length"]))]

    #Sort the data by type then by length
    data.sort(key=lambda x: (x["type"], x["length"]))

    types = [d["type"] for d in data]
    lengths = [d["length"] for d in data]
    distances = {key[1]: [d[key[0]] for d in data] for key in key_label}

    tick_labels = [f"{types[n]}_{lengths[n]}" for n in range(len(data))]

    #Get unique RNA types and their position
    unique_types, type_positions = np.unique(types, return_index=True)
    
    x = np.arange(len(data)) #Make locations for labels
    width = 0.005 #Bar width
    multiplier = 0

    fig, ax = plt.subplots(constrained_layout=True, figsize = (25, 12))

    #Plot groups
    for label, distance in distances.items(): 
        offset = width * multiplier
        ax.bar(x + offset, distance, width, label = label)
        multiplier += 1

    #Add vertical lines
    for position in type_positions[1:]: 
        ax.axvline(x=position - 0.5, linestyle = '--')

    #Set xticks and rotate
    ax.set_xticks(x + width * (multiplier-1)/2, tick_labels)
    ax.set_xticklabels(tick_labels, rotation = 90, fontsize = 16)

    ax.tick_params(axis='y', labelsize=18)
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}")) #Round y ticks to to decimal
    
    ax.legend(loc = 'upper left', bbox_to_anchor = (1.01, 1), fontsize = 16) #Set position of legend
    ax.set_ylabel("f", fontsize = 22)
    
    #Add horizontal grid lines behind bars
    ax.set_axisbelow(True)
    ax.grid(axis='y', linestyle = '--') 

    plt.savefig(output_file)
    return