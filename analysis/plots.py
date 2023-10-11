import matplotlib.pyplot as plt

### PLOT FUNCTIONS ###
def time_plot(): 
    return

def plot_times(data: dict, output_path):
    #NOTE - Copied from different file. Does not work as it is
    """
    Running time is plotted as a function of sequence length
    """
    x = data["length"] 
    y = data["time"]

    plt.scatter(x, y, color = 'black')
    plt.ylabel("Running time (s)")
    plt.xlabel("Sequence length")
    plt.grid(True)
    plt.savefig(os.path.join(output_path, "runningtimes.jpeg"))

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