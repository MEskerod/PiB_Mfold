import argparse, sys, math
import numpy as np
import pandas as pd
from Bio import SeqIO


def read_fasta(input) -> str:
    """
    Reads in a FASTA-file and returns the sequence
    If there is more than one sequence in the FASTA file it gives an error
    """
    sequences = []
    for record in SeqIO.parse(input, "fasta"): 
        sequences.append(record.seq)
    
    if len(sequences) > 1: 
        raise ValueError("FASTA file contains more than one sequence")

    return sequences[0]

def prepare_input(input: str) -> str: 
    """
    Strips new line in the end
    Makes sure that all letters are in upper case
    A value error is raised if any invalid letters are found in the sequence.  
    """
    allowed = ['A', 'C', 'G', 'U']
    input = input.strip().upper().replace('T', 'U')

    if any(letter not in allowed for letter in input): 
        raise ValueError("Invalid letters found in the sequence")
    
    return input

def read_parameters(file_loop, file_stacking): 
    """
    Read parameters from to .csv files:
    - One containing the loop parameters 
    - A second one containg the parameters for different types og base pairing 
    The .csv files are converted to Pandas tables and returned
    """
    #TODO - Update where loop parameters are used to handle loops_len>30
    loops = pd.read_csv(file_loop) 

    stacking = pd.read_csv(file_stacking, index_col=0)
    return loops, stacking

def loop_greater_10(loop_type, length, loop_parameters):
    """
    Calculates the energy parameters for loops with a size greater than 10 
    The parameter is calculated as described in 'Improved predictions of secondary structures for RNA'
    """
    max = {"IL":6, "BL":5, "HL":9}
    R = 0.001987 #In kcal/(mol*K)
    T = 310.15 #In K
    G_max = loop_parameters.at[max[loop_type],loop_type]

    G = G_max + 1.75*R*T*math.log(length/max[loop_type])

    return G

def find_E1(i, j, loop_parameters):
    """
    E1 is when Si and Sj basepair and gives one internal edge 
    This gives a hairpin loop 
    The function is able to handle loops of any size
    """
    size = j-i-1    
    energy = loop_parameters.at[size,"HL"]

    #NOTE - Below code alows loops of all sizes. Not used in example in article
    #if size <= 10: 
    #    energy = loop_parameters.at[size,"HL"]
    #else: 
    #    energy = loop_greater_10("HL", size, loop_parameters)

    return energy

def stacking(i, j, V, stacking_parameters, sequence): 
    """
    Find the energy parameter for basepairing of Sij and Si+1j-1, which results in basepair stacking
    If Si+1 and Sj+1 cannot basepair the energy is infinity
    Allows for Watson-Crick basepairs and wobble basepairs
    """
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

    prev_bp = sequence[i+1] + sequence[j-1]   
    
    if prev_bp in basepairs: 
        current_bp = sequence[i] + sequence[j]
        energy = stacking_parameters.at[current_bp, prev_bp] + V[i+1, j-1]
    else: 
        energy = float('inf')
    return energy

def bulge_loop(i, j, V, loop_parameters, stacking_parameters, sequence): 
    """
    Find the energy parameter of introducing a bulge loop on either "strand". 
    If the size of the bulge loop is one will basepairs on either side stack and stacking parameter is added
    Is able to handle loops of any size
    """
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']
    
    #Bulge on 3' end
    energy = float('inf')
    for jp in range(i+2,j-1):  
        bp = sequence[i+1]+sequence[jp]
        if bp in basepairs: 
            size = j-jp-1
            BL_energy = loop_parameters.at[size, "BL"] + V[i+1, jp]
            #If size = 1, stacking is retained
            if size == 1: 
                BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp]
            
            #NOTE - Code below allows loops of any size. Not used in example in article
            #if size <= 10:
            #   BL_energy = loop_parameters.at[size, "BL"] + V[i+1, jp]
               #If size = 1, stacking is retained
            #   if size == 1: 
            #       BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp] 
            #else: 
            #    BL_energy = loop_greater_10("BL", size, loop_parameters) + V[i+1, jp]
            

            if BL_energy < energy: 
                energy = BL_energy
    
    #Bulge on 5' end
    for ip in range(i+2,j-1):  
        bp = sequence[ip]+sequence[j-1]
        if bp in basepairs: 
            size = ip-i-1
            BL_energy = loop_parameters.at[size, "BL"] + V[ip, j-1]
            #If size = 1, stacking is retained
            if size == 1: 
                BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp]
            
            #NOTE - Code below allows loops of any size. Not used in example in article
            #if size <= 10:
            #   BL_energy = loop_parameters.at[size, "BL"] + V[ip, j-1]
               #If size = 1, stacking is retained
            #   if size == 1: 
            #       BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp] 
            #else: 
            #    BL_energy = loop_greater_10("BL", size, loop_parameters) + V[ip, j-1]

            if BL_energy < energy: 
                energy = BL_energy
    return energy

def interior_loop(i, j, V, loop_parameters, sequence): 
    """
    Find the energy parameter of adding a interior loop. 
    Penalty is added for asymmetric loops.
    Is able to handle loops of any size
    """
    basepairs = ['CG', 'GC', 'GU', 'UG']
    penalty_pairs = ['AU', 'UA']
    
    energy = float('inf')
    
    bp = sequence[i]+sequence[j]

    for ip in range(i+2, j-2): 
        for jp in range(ip+3, j-1): 
            bp_prime = sequence[ip] + sequence[jp]
            if bp_prime in basepairs+penalty_pairs:
                size = (ip-i-1)+(j-jp-1)
                IL_energy = loop_parameters.at[size, "IL"] + V[ip, jp]

                #NOTE - Below code for loops of any length. Not used in example in article
                #if size <= 10: 
                #    IL_energy = loop_parameters.at[size, "IL"] + V[ip, jp]
                #else: 
                #    IL_energy = loop_greater_10("IL", size, loop_parameters) + + V[ip, jp]
                
                
                #NOTE - Penalty for asymmetric loops described in artice, but not used in example
                #Add penalty to energy if loop is asymmetric
                #if (ip-i-1) != (j-jp-1): 
                #    f = [0.4, 0.3, 0.2, 0.1]
                #    N1 = (ip-i-1)
                #    N2 =(j-jp-1)
                #    N = abs(N1-N2)
                #    M = min(4, N1, N2)-1
                #    penalty = min(3, N*f[M]) 
                #    IL_energy += penalty
                
                #Check if energy is smaller than current min
                if IL_energy < energy: 
                    energy = IL_energy
    return energy

def find_E2(i, j, V, parameters, sequence): 
    """
    E2 is when Si and Sj contains two internal edged 
    Contains an edge between Si and Sj and an edge between Si' and Sj' 
    i<i'<j'<j
    Can be stacking, bulge loops or internal loops
    Returns the minimum of the 3 options. 
    """
    loop_parameters =  parameters[0]
    stacking_parameters = parameters[1]

    energy = min(stacking(i, j, V, stacking_parameters, sequence), bulge_loop(i, j, V, loop_parameters, stacking_parameters, sequence), interior_loop(i, j, V, loop_parameters, sequence))
    return energy

def find_E3(i, j, W): 
    """
    E3 contains more than two internal edges
    Gives a bifurcating loop
    The energy is the energy of the substructures 
    i+1<i'<j-2
    """
    energy = float('inf')
    for ip in range(i+2, j-2):  #NOTE - Dobbelt check that no base pairing condition is needed
        loop_energy = W[i+1, ip] + W[ip+1, j-1]
        if loop_energy < energy: 
            energy = loop_energy
    return energy

def penta_nucleotides(W, V, sequence, parameters):
    """
    Initiates the V and W matrices. 
    The shortest possible subsequences are of length 5 and can only form hairpin loops of size 3 if i and j basepair
    """
    N = len(sequence)
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

    loop_parameters = parameters[0]
    
    for i in range(0, N-4): 
        j = i+4
        if sequence[i]+sequence[j] not in basepairs:
            V[i,j] = W[i,j ]= float('inf')
        else: 
            V[i,j] = W[i,j] = loop_parameters.at[3, "HL"] 

def compute_V(i, j, W, V, sequence, parameters): 
    """
    Computes the minimization over E1, E2 and E3, which will give the value at V[i,j]
    """
    
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

    if sequence[i] + sequence[j] in basepairs:
        v = round(min(find_E1(i, j, parameters[0]), find_E2(i, j, V, parameters, sequence), find_E3(i, j, W)),5)
    else: 
        v = float('inf')

    V[i, j] = v

def find_E4(i, j, W): 
    """
    i and j both base pair, but not with each other. 
    It find the minimum of combinations of to possible subsequences containing i and j
    """
    energy = float('inf')
    for ip in range(i+1, j-1): 
        E = W[i, ip] + W[ip+1, j]
        if E < energy: 
            energy = E
    return energy

def compute_W(i, j, W, V):
    """
    Computes the minimization over possibilities for W, which will give the value for W[i,j]
     Possibilities are: 
    - i or j in a structure (W[i+1, j] or W[i, j-1])
    - i and j basepair with each other (V[i,j])
    - i and j both base pair but not with each other (E4)
    """
    w = min(W[i+1,j], 
            W[i,j-1], 
            V[i,j], 
            find_E4(i, j, W))
    W[i,j] = w

def fold_rna(sequence, parameters): 
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
    penta_nucleotides(W, V, sequence, parameters) 

    for l in range(5, N): #Computes the best score for all subsequences that are longer than 5 nucleotides
        for i in range(0, N-5): 
            j = i+l
            if j < N: 
                compute_V(i, j, W, V, sequence, parameters)
                compute_W(i, j, W, V)
    
    return W, V

def find_optimal(W) -> float: 
    """
    Find the final energy of the folded RNA
    """
    return W[0, -1]

def backtrack(W, V): 
    """
    Backtracks trough the W, V matrices to find the final fold
    """
    #TODO - Finish backtracking
    #NOTE - Is desired output dot bracket?
    dotbracket = ""




    return dotbracket

def main() -> None: 
    """
    Running the program
    Takes following arguments: 
    input - Can be a file or written in the command line. -i followed by sequence or -f followed by file name. File must be a FASTA file containing one sequence
    ouput - Can be a file specified by using the flag -o or a the default stdout
    """
    #Setting up the option parsing using the argparse module
    argparser = argparse.ArgumentParser(
        description="" )
    #Adding arguments
    #TODO - Add description/help for command line options
    #Input can either be provided in a file or in stdin
    argparser.add_argument('-i', '--input') 
    argparser.add_argument('-f', '--file', type=argparse.FileType('r'))
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('-o', '--outfile', metavar='output', type=argparse.FileType('w'), default=sys.stdout)

    args = argparser.parse_args()

    if args.input: 
        sequence = prepare_input(args.input)

    if args.file: 
        sequence = prepare_input(read_fasta(args.file))

    
    #TODO - Finish main
    parameters = read_parameters("loop_parameters.csv", "pairing_parameters.csv")
    W, V = fold_rna(sequence, parameters)
    energy = find_optimal(W)
    #fold = backtrack(W, V, sequence)

    
    #NOTE - Just for "testing"
    print("V:\n", V)
    print("W:\n", W)
    print("Energy of optimal fold: ", energy)
    

    return #TODO - Figure out what the desired output is. Fold? Fold + energy?

if __name__ == '__main__': 
    main()