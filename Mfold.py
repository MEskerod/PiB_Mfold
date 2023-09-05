import argparse, sys, fileinput
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
    #TODO - Error if len(sequences) > 1

    return sequences[0]

def prepare_input(input: list[str]) -> str: 
    """
    Removes name of sequence if input is fasta 
    Strips any whitespace
    Makes sure that all letters are in upper case 
    """
    input = "".join([line.strip().upper() for line in input if line[0].isalpha()])
    return input

def read_parameters(file_loop, file_stacking): 
    """
    """
    loops = pd.read_csv(file_loop)

    stacking = pd.read_csv(file_stacking, index_col=0)
    return loops, stacking


def find_E1(i, j, parameters):
    """
    E1 is when Si and Sj basepair and gives one internal edge 
    This gives a hairpin loop 
    """
    loop_parameters = parameters[0]
    energy = loop_parameters.at[(j-i-1),"HL"]
    return energy

def stacking(i, j, V, stacking_parameters, sequence): 
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

    prev_bp = sequence[i+1] + sequence[j-1] 
    
    if prev_bp in basepairs: 
        current_bp = sequence[i] + sequence[j]
        energy = stacking_parameters.at[current_bp, prev_bp] + V[i+1, j-1]
    else: 
        energy = float('inf')
    return energy

def bulge_loop(i, j, V, loop_parameters, stacking_parameters, sequence): 
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']
    
    energy = float('inf')
    for jp in range(i+2,j-1):  
        bp = sequence[i+1]+sequence[jp]
        if bp in basepairs: 
            size = j-jp-1
            BL_energy = loop_parameters.at[size, "BL"] + V[i+1, jp]
            #If size = 1, stacking is retained
            if size == 1: 
                BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp]
            if BL_energy < energy: 
                energy = BL_energy
    
    for ip in range(i+2,j-1):  
        bp = sequence[ip]+sequence[j-1]
        if bp in basepairs: 
            size = ip-i-1
            BL_energy = loop_parameters.at[size, "BL"] + V[ip, j-1]
            #If size = 1, stacking is retained
            if size == 1: 
                BL_energy += stacking_parameters.at[(sequence[i]+sequence[j]), bp]
            if BL_energy < energy: 
                energy = BL_energy
    return energy

def interior_loop(i, j, V, loop_parameters, sequence): 
    """
    """
    basepairs = ['CG', 'GC', 'GU', 'UG']
    penalty_pairs = ['AU', 'UA']
    
    energy = float('inf')
    
    bp = sequence[i]+sequence[j]

    for ip in range(i+2, j-2): 
        for jp in range(ip+3, j-1): 
            bp_prime = sequence[ip] + sequence[jp]
            if bp_prime in basepairs+penalty_pairs:
                #If both are in basepairs
                if all(pair in basepairs for pair in [bp, bp_prime]):
                    size = (ip-i-1)+(j-jp-1)
                    IL_energy = loop_parameters.at[size, "IL"] + V[ip, jp]
                #If one of the pairs is AU penalty is added
                elif (bp in basepairs and bp_prime in penalty_pairs) or (bp in penalty_pairs and bp_prime in basepairs): 
                    size = (ip-i-1)+(j-jp-1)
                    IL_energy = loop_parameters.at[size, "IL"] + V[ip, jp] #+ 0.9 #NOTE - Penalty for AU base pairs described in article but not used in example
                #If both are AU or GU doubbel penalty is added
                elif all(pair in penalty_pairs for pair in [bp, bp_prime]):
                    size = (ip-i-1)+(j-jp-1)
                    IL_energy = loop_parameters.at[size, "IL"] + V[ip, jp] #+ 1.8
                
                #NOTE - Penalty described in artice, but not used in example
                #Add penalty to energy if loop is asymmetric
                #if (ip-i-1) != (j-jp-1): 
                    #f = [0.7, 0.6, 0.4, 0.2, 0.1]
                    #N1 = (ip-i-1)
                    #N2 =(j-jp-1)
                    #N = abs(N1-N2)
                    #M = min(5, N1, N2)-1
                    #penalty = min(6, N*f[M]) 
                    #IL_energy += penalty
                
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
    Computes the minimization over E1, E2 and E3
    """
    
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']

    if sequence[i] + sequence[j] in basepairs:
        v = round(min(find_E1(i, j, parameters), find_E2(i, j, V, parameters, sequence), find_E3(i, j, W)),5)
    else: 
        v = float('inf')

    V[i, j] = v

def find_E4(i, j, W): 
    """
    i and j both base pair, but not with each other
    """
    energy = float('inf')
    for ip in range(i+1, j-1): 
        E = W[i, ip] + W[ip+1, j]
        if E < energy: 
            energy = E
    return energy

def compute_W(i, j, W, V):
    """
    Computes the minimization over possibilities for W
    """
    w = min(W[i+1,j], 
            W[i,j-1], 
            V[i,j], 
            find_E4(i, j, W))
    W[i,j] = w

def fold_rna(sequence, parameters): 
    """
    Fills out the W and V matrices to find the fold that gives the minimum free energy
    Folows Mfold
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
    dotbracket = ""




    return dotbracket

def main() -> None: 
    """
    Running the program
    Takes following arguments: 
    input - Can be a file or written in the command line. -i followed by sequence or -f followed by file name. File must be a FASTA file
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

    parameters = read_parameters("loop_parameters.csv", "pairing_parameters.csv")
    
    #TODO - Finish main

    #print(sequence)

    W, V = fold_rna(sequence, parameters)

    energy = find_optimal(W)

    #fold = backtrack(W, V, sequence)

    print("V:\n", V)
    print("W:\n", W)
    print("Energy of optimal fold: ", energy)


    

    return

if __name__ == '__main__': 
    main()