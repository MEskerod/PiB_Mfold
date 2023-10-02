import math
import numpy as np
import pandas as pd

### FUNCTIONS USED BY SEVERAL VERSIONS ### 
def loop_greater_10(max_lengths: dict, loop_parameters: pd.array):
    """
    Calculates the energy parameters for loops with a size greater than 10 
    The parameter is calculated as described in 'Improved predictions of secondary structures for RNA'

    The max_lengths is a dictionary contain the length of the longest sequence of each type that's exeperimentally determined: 
    max_lenghts = {"IL":?, "BL":?, "HL":?}

    The function returns a function that can calculate loop energies given length, type
    """
    def loop_function(loop_type, length):
        R = 0.001987 #In kcal/(mol*K)
        T = 310.15 #In K
        G_max = loop_parameters.at[max_lengths[loop_type],loop_type]

        G = G_max + 1.75*R*T*math.log(length/max_lengths[loop_type])

        return G
    
    return loop_function

def asymmetric_penalty(f, penalty_max): 
    """
    f has to be a list. In articles writen as f(1), f(2) and so forth

    The asymmetry function is used to calculate a penalty to internal loops that are asymmetric. 
    This penalty does not exists in the orginial paper, but is added later

    This functions returns a function that uses the given parameters to calculate the penalty for asymmetric loops of given size
    """

    M_max = len(f)

    def asymmetry_function(i, ip, j, jp):
        N1 = (ip-i-1)
        N2 =(j-jp-1)
        N = abs(N1-N2)
        M = min(M_max, N1, N2)-1
        penalty = min(penalty_max, N*f[M]) 
    
        return penalty
    
    return asymmetry_function



### LOOP ENERGIES ###
def stacking(i, j, V, stacking_parameters, sequence): 
    """
    Find the energy parameter for basepairing of Sij and Si+1j-1, which results in basepair stacking
    If Si+1 and Sj+1 cannot basepair the energy is infinity
    Allows for Watson-Crick basepairs and wobble basepairs
    """
    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}

    prev_bp = sequence[i+1] + sequence[j-1]   
    
    if prev_bp in basepairs: 
        current_bp = sequence[i] + sequence[j]
        energy = round(stacking_parameters.at[current_bp, prev_bp] + V[i+1, j-1], 5)
    else: 
        energy = float('inf')
    return energy

def bulge_loop_3end(i, j, V, loop_parameters, stacking_parameters, sequence, stacking: bool): 
    """
    Find the energy parameter of introducing a bulge loop on the 3' end. 
    If the size of the bulge loop is 1 and stacking=True, stacking on each side of the loop is retained and stacking parameter is added
    Is able to handle loops of any size
    """
    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}
    energy = float('inf')
    ij = None

    #Bulge on 3' end
    for jp in range(i+2,j-1):  
        bp = sequence[i+1]+sequence[jp]
        if bp in basepairs: 
            size = j-jp-1
            if size <= 10:
               BL_energy = loop_parameters.at[size, "BL"] + V[i+1, jp]
               if size == 1 and stacking: 
                   BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp]
            else: 
                BL_energy = loop_greater_10("BL", size, loop_parameters) + V[i+1, jp]
            
            if BL_energy < energy: 
                energy = round(BL_energy, 5)
                ij = jp
    
    return energy, ij

def bulge_loop_5end(i, j, V, loop_parameters, stacking_parameters, sequence, stacking: bool):
    """
    Find the energy parameter of introducing a bulge loop on the 5' end. 
    If the size of the bulge loop is 1 and stacking=True, stacking on each side of the loop is retained and stacking parameter is added
    Is able to handle loops of any size
    """
    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}
    energy = float('inf')
    ij = None
    
    #Bulge on 5' end
    for ip in range(i+2,j-1):  
        bp = sequence[ip]+sequence[j-1]
        if bp in basepairs: 
            size = ip-i-1
            if size <= 10:
                BL_energy = loop_parameters.at[size, "BL"] + V[ip, j-1]
                if size == 1 and stacking: 
                    BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp] 
            else: 
                BL_energy = loop_greater_10("BL", size, loop_parameters) + V[ip, j-1]

            if BL_energy < energy: 
                energy = round(BL_energy, 5)
                ij = ip

    return energy, ij

def interior_loop(i, j, V, loop_parameters, sequence, closing_penalty: bool, asymmetry: bool): 
    """
    Find the energy parameter of adding a interior loop. 
    Is able to handle loops of any size
    A penalty is added for asymmetric loops. If the penalty should not be added interior_loop should be called with asymmetry = False
    A penalty is added for interior loops closed by AU and GU base pairs
    """
    basepairs = {'CG', 'GC', 'GU', 'UG', 'AU', 'UA'}
    
    energy = float('inf')
    ij = None

    for ip in range(i+2, j-2): 
        for jp in range(ip+3, j-1): 
            bp_prime = sequence[ip] + sequence[jp]
            if bp_prime in basepairs:
                size = (ip-i-1)+(j-jp-1)

                IL_energy = (loop_parameters.at[size, "IL"] + V[ip, jp]) if size <= 10 else (loop_greater_10("IL", size, loop_parameters) + V[ip, jp])
                
                #Add penalty to energy if loop is asymmetric
                if asymmetry and ((ip-i-1) != (j-jp-1)): 
                    IL_energy += asymmetric_penalty(i, ip, j, jp)

                #Add penalty if closing base pairs are AU og GU base pairs
                if closing_penalty: 
                    bp = sequence[i] + sequence[j]
                    if bp in ['GU', 'UG', 'AU', 'UA']: 
                        IL_energy += 0.9
                    if bp_prime in ['GU', 'UG', 'AU', 'UA']:
                        IL_energy += 0.9 
                
                #Check if energy is smaller than current min
                if IL_energy < energy: 
                    energy = round(IL_energy, 5)
                    ij = (ip, jp)
    
    return energy, ij

def find_E1(i, j, loop_parameters):
    """
    E1 is when Si and Sj basepair and gives one internal edge 
    This gives a hairpin loop 
    The function is able to handle loops of any size
    """
    size = j-i-1    

    energy = loop_parameters.at[size,"HL"] if size <= 10 else loop_greater_10("HL", size, loop_parameters)

    return round(energy, 5)

def find_E2(i, j, V, parameters, sequence, stacking: True, closing_penalty: True, asymmetry: True): 
    """
    E2 is when Si and Sj contains two internal edged 
    Contains an edge between Si and Sj and an edge between Si' and Sj' 
    i<i'<j'<j
    Can be stacking, bulge loops or internal loops
    Returns the minimum of the 3 options. 
    """
    loop_parameters =  parameters[0]
    stacking_parameters = parameters[1]

    energy = min(stacking(i, j, V, stacking_parameters, sequence), 
                 bulge_loop_3end(i, j, V, loop_parameters, stacking_parameters, sequence, stacking)[0], 
                 bulge_loop_5end(i, j, V, loop_parameters, stacking_parameters, sequence, stacking)[0], 
                 interior_loop(i, j, V, loop_parameters, sequence, closing_penalty, asymmetry)[0])
    return energy

def find_E3(i, j, W): 
    #TODO - Check if this one has to be changed to include other options for bifurcating loops
    """
    E3 contains more than two internal edges
    Gives a bifurcating loop
    The energy is the energy of the substructures 
    i+1<i'<j-2
    """
    energy = float('inf')
    ij = None

    for ip in range(i+2, j-2):  
        loop_energy = W[i+1, ip] + W[ip+1, j-1]
        if loop_energy < energy: 
            energy = round(loop_energy, 5)
            ij = (ip, ip+1)
    return energy, ij

def find_E4(i, j, W): 
    """
    i and j both base pair, but not with each other. 
    It find the minimum of combinations of to possible subsequences containing i and j
    """
    energy = float('inf')
    ij = None

    for ip in range(i+1, j-1): 
        subsequence_energy = W[i, ip] + W[ip+1, j]
        
        if subsequence_energy < energy: 
            energy = round(subsequence_energy, 5)
            ij = (ip, ip+1)

    return energy, ij

### FILL V AND W ###
def penta_nucleotides(W, V, sequence, loop_parameters):
    """
    Initiates the V and W matrices. 
    The shortest possible subsequences are of length 5 and can only form hairpin loops of size 3 if i and j basepair
    """
    N = len(sequence)
    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}
    
    for i in range(0, N-4): 
        j = i+4
        bp = sequence[i]+sequence[j]
        if bp not in basepairs:
            V[i,j] = W[i,j ]= float('inf')
        else: 
            V[i,j] = W[i,j] = loop_parameters.at[3, "HL"] 