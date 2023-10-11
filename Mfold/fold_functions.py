import math
import numpy as np
import pandas as pd

#TODO - Read and update all descriptons to match the content!
#TODO - Add types for all parameters

### FUNCTIONS USED BY SEVERAL VERSIONS ### 

def make_asymmetric_penalty(f, penalty_max): 
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

def loop_greater_10(loop_type, length, loop_parameters: pd.array):
    """
    Calculates the energy parameters for loops with a size greater than 10 
    The parameter is calculated as described in 'Improved predictions of secondary structures for RNA'

    The function returns the energi parameter for a loop of a given type and size. 
    """
    R = 0.001987 #In kcal/(mol*K)
    T = 310.15 #In K
    G_max = loop_parameters.at[10, loop_type]

    G = G_max + 1.75*R*T*math.log(length/10)

    return G

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

def interior_loop(i, j, V, loop_parameters, sequence, asymmetric_penalty_function, closing_penalty: bool, asymmetry: bool): 
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
                    IL_energy += asymmetric_penalty_function(i, ip, j, jp)

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

def find_E2(i, j, V, parameters, sequence, asymmetric_penalty_function, bulge_stacking: bool, closing_penalty: bool, asymmetry_penalty: bool): 
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
                 bulge_loop_3end(i, j, V, loop_parameters, stacking_parameters, sequence, bulge_stacking)[0], 
                 bulge_loop_5end(i, j, V, loop_parameters, stacking_parameters, sequence, bulge_stacking)[0], 
                 interior_loop(i, j, V, loop_parameters, sequence, asymmetric_penalty_function, closing_penalty, asymmetry_penalty)[0])
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

def compute_V(i, j, W, V, sequence, parameters, asymmetric_penalty_function, bulge_stacking: bool, closing_penalty: bool, asymmetry_penalty: bool): 
    """
    Computes the minimization over E1, E2 and E3, which will give the value at V[i,j]
    """
    
    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}

    if sequence[i] + sequence[j] in basepairs:
        v = min(find_E1(i, j, parameters[0]), 
                find_E2(i, j, V, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty), 
                find_E3(i, j, W)[0])

    else: 
        v = float('inf')

    V[i, j] = v

def compute_W(i, j, W, V):
    """
    Computes the minimization over possibilities for W, which will give the value for W[i,j]
     Possibilities are: 
    - i or j in a structure (W[i+1, j] or W[i, j-1])
    - i and j basepair with each other (V[i,j])
    - i and j both base pair but not with each other (E4)
    """
    w = min(W[i+1,j], W[i,j-1], V[i,j], find_E4(i, j, W)[0])

    W[i,j] = w


def fold_rna(sequence, parameters, asymmetric_penalty_function, bulge_stacking: bool, closing_penalty: bool, asymmetry_penalty: bool): 
    """
    Fills out the W and V matrices to find the fold that gives the minimum free energy
    Follows Mfold as desribed by M. Zuker

    The V matrix contains the minimum free energy for the subsequences i and j, if i and j has to form a pair. 
    Is Si and Sj are not able to basepair the energy will be infinity (not possible)

    The W matrix contains the minimum free energy for the subsequences i and j, but i and j does not have to basepair. 
    """
    N = len(sequence)
    W, V = np.full([N, N], float('inf')), np.full([N, N], float('inf'))


    #Fills out the table with all the penta nucleotide.
    #Penta nucleotides are the base cases. If shorter they cannot be folded
    penta_nucleotides(W, V, sequence, parameters[0]) 

    for l in range(5, N): #Computes the best score for all subsequences that are longer than 5 nucleotides
        for i in range(0, N-5): 
            j = i+l
            if j < N: 
                compute_V(i, j, W, V, sequence, parameters, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty) 
                compute_W(i, j, W, V)

    return W, V

def find_optimal(W) -> float: 
    """
    Find the final energy of the folded RNA
    """
    return W[0, -1]

### BACTRACKING ### 
def trace_V(i, j, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking: bool, closing_penalty: bool, asymmetry_penalty: bool): 
    """
    """
    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}
    
    loop_parameters, stacking_parameters = parameters[0], parameters[1]
    
    if V[i,j] == find_E1(i, j, loop_parameters): 
        dotbracket[i], dotbracket[j] = '(', ')'
        for n in range(i+1, j): 
            dotbracket[n] = '.'
    
    elif V[i,j] == stacking(i, j, V, stacking_parameters, sequence): 
        dotbracket[i], dotbracket[j] = '(', ')'
        trace_V(i+1, j-1, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)
    
    elif V[i,j] == bulge_loop_3end(i, j, V, loop_parameters, stacking_parameters, sequence, bulge_stacking)[0]: 
        jp = bulge_loop_3end(i, j, V, loop_parameters, stacking_parameters, sequence, bulge_stacking)[1]
        dotbracket[i], dotbracket[j] = '(', ')'
        for n in range(jp, j): 
            dotbracket[n] = '.'
        trace_V(i+1, jp, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)
    
    elif V[i,j] == bulge_loop_5end(i, j, V, loop_parameters, stacking_parameters, sequence, bulge_stacking)[0]: 
        ip = bulge_loop_5end(i, j, V, loop_parameters, stacking_parameters, sequence, bulge_stacking)[1]
        dotbracket[i], dotbracket[j] = '(', ')'
        for n in range(i+1, ip): 
            dotbracket[n] = '.'
        trace_V(ip, j-1, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)
    
    elif V[i,j] == interior_loop(i, j, V, loop_parameters, sequence, asymmetric_penalty_function, closing_penalty, asymmetry_penalty)[0]:
        ij = interior_loop(i, j, V, loop_parameters, sequence, asymmetric_penalty_function, closing_penalty, asymmetry_penalty)[1]
        dotbracket[i], dotbracket[j] = '(', ')' 
        for n in range(i+1, ij[0]): 
            dotbracket[n] = '.'
        for n in range(ij[1]+1, j): 
            dotbracket[n] = '.'
        trace_V(ij[0], ij[1], W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)
    
    elif V[i, j] == find_E3(i, j, W)[0]: 
        ij = find_E3(i, j, W)[1]
        dotbracket[i], dotbracket[j] = '(', ')' 
        trace_W(i+1, ij[0], W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty), trace_W(ij[1], j-1, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)

def trace_W(i, j, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking: bool, closing_penalty: bool, asymmetry_penalty: bool): 
    """
    """
    if W[i,j] == W[i+1, j]: 
        dotbracket[i] = '.'
        trace_W(i+1, j, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)

    elif W[i,j] == W[i, j-1]: 
        dotbracket[j] = '.'
        trace_W(i, j-1, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)

    elif W[i, j] == V[i, j]: 
        trace_V(i, j, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)

    elif W[i,j] == find_E4(i, j, W)[0]: 
        ij = find_E4(i,j,W)[1] 
        trace_W(i, ij[0], W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty), trace_W(ij[1], j, W, V, dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)



def backtrack(W, V, parameters, sequence, asymmetric_penalty_function, bulge_stacking: bool, closing_penalty: bool, asymmetry_penalty: bool): 
    """
    Backtracks trough the W, V matrices to find the final fold
    """
    dotbracket =  ['?' for x in range(W.shape[0])]
    
    j = W.shape[0]-1
    i = 0
    
    trace_W(i, j, W, V,dotbracket, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)

    return "".join(dotbracket)