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

def read_parameters(file): 
    """
    """
    loops = pd.read_csv(file)
    return loops

def base_pairs() -> list: 
    """
    """
    #NOTE - Delete if not usefull
    return [('A', 'U'), ('U', 'A'), ('G', 'C'), ('C', 'G')] 

def find_E1(i, j, parameters, sequence): 
    """
    E1 is when Si and Sj basepair and gives one internal edge 
    This gives a hairpin loop 
    """
    return

def find_E2(i, j, V, parameters, sequence): 
    """
    E2 is when Si and Sj contains two internal edged 
    Contains an edge between Si and Sj and an edge between Si' and Sj' 
    i<i'<j'<j
    """
    return

def find_E3(i, j, W): 
    """
    E3 contains more than two internal edges
    Gives a bifurcating loop
    The energy is the energy of the substructures 
    i+1<i'<j-2
    """
    return

def penta_nucleotides(W, V, sequence, parameters):
    """
    """
    N = len(sequence)
    basepairs = base_pairs()
    
    for i in range(0, N-4): 
        j = i+4
        if (sequence[i], sequence[j]) not in basepairs:
            V[i,j] = float('inf')
        else: 
            V[i,j] = parameters.at[3, "HL"]

def compute_V(i, j, W, V, sequence, parameters): 
    """
    Computes the minimization over E1, E2 and E3
    """
    basepairs = base_pairs
    if (sequence[i], sequence[j]) not in basepairs: 
        v = float('inf')
    else:
        v = min(find_E1(i, j, parameters, sequence), find_E2(i, j, V, parameters, sequence), find_E3(i, j, W))
    V[i, j] = v

def find_E4(i, j, W): 
    """
    i and j both base pair, but not with each other
    """
    return

def compute_W(i, j, W, V):
    """
    Computes the minimization over possibilities for W
    """
    w = min(W[i+1,j], W[i,j-1], V[i,j], find_E4(i, j, W))
    W[i,j] = w

def fold_rna(sequence, parameters): 
    """
    Fills out the W and V matrices to find the fold that gives the minimum free energy
    Folows Mfold
    """
    N = len(sequence)
    W, V = np.zeros([N, N]), np.zeros([N, N])

    #TODO - update/check for loop
    #for i in range(0, N)

    penta_nucleotides(W, V, sequence, parameters)

    #for l in range(5, N): #We starts by computing all pentanucleotides and continues with larger sequences
    #    for i in range(0, N-5): 
    #        j = i+l
    #        if j < N: 
                #compute_W(i, j, W, V)
                #compute_V(i, j, W, V, sequence, parameters) 
    
    return W, V

def find_optimal(W) -> float: 
    """
    Find the final energy of the folded RNA
    """
    return W[-1, -1]

def backtrack(): 
    """
    Backtracks trough the W, V matrices to find the final fold
    """
    return

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
    #TODO - Set up arguments
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

    loop_parameters = read_parameters("loop_parameters.csv")
    
    #print(sequence)

    W, V = fold_rna(sequence, loop_parameters)

    #energy = find_optimal(W)

    #fold = backtrack(W, V, sequence)

    print(V)

    

    return

if __name__ == '__main__': 
    main()