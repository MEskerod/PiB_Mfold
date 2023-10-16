import matplotlib.pyplot as plt
import os

### PLOT FUNCTIONS ###
def plot_synthetic_real_times(lengths, synthetic_time, real_time, output_path):
    """
    Running time is plotted as a function of sequence length
    """
    plt.scatter(lengths, real_time, label = "Real RNA sequences", s=20, edgecolors='blue', facecolors = 'none', linewidths=1)
    plt.plot(lengths, real_time, color = 'blue', linestyle = '--', linewidth = 0.8)
    plt.scatter(lengths, synthetic_time, label = "Synthetic sequences", s=20, edgecolors= 'red', facecolors = 'none', linewidths=1)
    plt.plot(lengths, synthetic_time, color = 'red', linestyle = '--', linewidth = 0.8)
    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.legend()
    plt.grid(True, linestyle = '--')
    plt.savefig(os.path.join(output_path, "synthetic_vs_real_times.jpeg"))

def plot_Nussinov(): 
    """
    """
    return

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