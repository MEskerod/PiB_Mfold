import argparse, sys, time
from io import TextIOWrapper
from Bio import SeqIO
import numpy as np

def read_fasta(input: str) -> str:
    """
    Reads in a FASTA-file and returns the sequence
    If there is more than one sequence in the FASTA file it gives an error
    """
    records = list(SeqIO.parse(input, 'fasta'))
    
    if len(records) > 1: 
        raise ValueError("FASTA file contains more than one sequence")

    return str(records[0].seq), records[0].id

def prepare_input(input: str) -> str: 
    """
    Strips new line in the end
    Makes sure that all letters are in upper case
    A value error is raised if any invalid letters are found in the sequence.  
    """
    allowed = set(['A', 'C', 'G', 'U', 'N'])
    input = input.strip().upper().replace('T', 'U')

    if any(letter not in allowed for letter in input): 
        raise ValueError("Invalid letters found in the sequence")
    
    return input

def db_to_file(sequence: str, db: str, filename: str, name: str) -> None:
    """
    Writes sequence and structure to a dbn file. 
    Header lines are denoted by #
    """ 
    with open(filename, 'w') as f: 
        f.write(f"#Name: {name}\n")
        f.write(f"#Length: {len(sequence)}\n")
        f.write(sequence + "\n")
        f.write(db + "\n")

def pairing_score(i: int, j: int, S: np.array, sequence: str): 
    """
    Returns the score of pairing i and j. 
    If i and j can form a basepair the score is 1 + S[i+1, j-1] 
    Otherwise it is negative infinity (not possible structure)
    """
    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}

    if sequence[i]+sequence[j] in basepairs: 
        score = 1 + S[i+1, j-1]
    else: 
        score = float('-inf')
    
    return score

def bifurcating_score(i: int, j: int, S: np.array) -> tuple:
    """
    Tries all the possible bifurcating loops and returns the score and k that gives the maximum energy
    """
    score = float('-inf')
    end_k = 0

    for k in range(i+1, j-1): 
        sub_score = S[i, k] + S[k+1, j]
        if sub_score > score: 
            score = sub_score
            end_k = k
    return score, end_k

def fill_S(sequence: str) -> np.array:
    """
    Fills out the S matrix
    """

    N = len(sequence)
    S = np.zeros([N, N])

    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}

    
    for l in range(4, N): #Computes the best score for all subsequences that are 5 nucleotides or longer
        for i in range(0, N-4): 
            j = i+l
            if j < N:
                score = max(S[i+1, j],
                            S[i, j-1],
                            pairing_score(i, j, S, sequence),
                            bifurcating_score(i, j, S)[0])
                S[i,j] = score
                
    return S

def find_optimal(S: np.array) -> float: 
    """
    Find the final energy of the folded RNA
    """
    return S[0, -1]

def backtrack(i: int, j: int, S: np.array, dotbracket: list, sequence: str) -> None: 
    """
    Backtracks trough the S matrix, to find the structure that gives the maximum energy
    """
    if j-i-1 <= 3: 
        dotbracket[i], dotbracket[j] = '(', ')'
        for n in range(i+1, j): 
            dotbracket[n] = '.'

    elif S[i, j] == S[i+1, j]: 
        dotbracket[i] = '.'
        backtrack(i+1, j, S, dotbracket, sequence)

    elif S[i, j] == S[i, j-1]: 
        dotbracket[j] = '.'
        backtrack(i, j-1, S, dotbracket, sequence)
    elif S[i, j] == pairing_score(i, j, S, sequence): 

        dotbracket[i], dotbracket[j] = '(', ')'
        backtrack(i+1, j-1, S, dotbracket, sequence)

    elif S[i, j] == bifurcating_score(i, j, S)[0]:
        k = bifurcating_score(i, j, S)[1]
        backtrack(i, k, S, dotbracket, sequence), backtrack(k+1, j, S, dotbracket, sequence)

    

def fold_RNA(S: np.array, sequence: str) -> str: 
    """
    Finds the optimal structure of foldning the sequence and returns the dot bracket notation
    """
    dotbracket =  ['?' for x in range(S.shape[0])]
    
    j = S.shape[0]-1
    i = 0

    backtrack(i, j, S, dotbracket, sequence)

    return "".join(dotbracket)


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
    #Input can either be provided in a file or in stdin
    argparser.add_argument('-i', '--input', metavar='', help = 'input provided at sequence in command line') 
    argparser.add_argument('-f', '--file', type=argparse.FileType('r'), metavar='', help = 'input provided as fasta file')
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('-o', '--outfile', metavar='', default=sys.stdout, help='name for output file without extension')

    args = argparser.parse_args()

    if args.input: 
        sequence = prepare_input(args.input)
        name = "User inputted sequence"

    if args.file: 
        sequence, name = read_fasta(args.file)
        sequence = prepare_input(sequence)
    
    if not sequence: 
        raise ValueError("No valid input sequence provided.")

    print(f"Fold {name}")
    start_time = time.time()

    S = fill_S(sequence)
    energy = find_optimal(S)
    fold = fold_RNA(S, sequence)

    #Write to outfile
    if args.outfile == sys.stdout: 
        args.outfile.write(fold + "\n\n")
    
    else: 
        outfile = args.outfile + ".dbn"
        db_to_file(sequence, fold, outfile, name)
    
    print("Finished in {} second".format(round(time.time() - start_time, 2)))
    print(f"Energy of optimal fold is {energy}\n")

   
if __name__ == '__main__': 
    main()