import ast, argparse, sys, time, math, os
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import TextIOWrapper

####### HELP FUNCTIONS #######
def read_fasta(input: str) -> str:
    """
    Reads in a FASTA-file and returns the sequence
    If there is more than one sequence in the FASTA file an error is raised 

    Args: 
        input: path to input sequence
    """
    records = list(SeqIO.parse(input, 'fasta'))
    
    if len(records) > 1: 
        raise ValueError("FASTA file contains more than one sequence")

    return str(records[0].seq), records[0].id

def prepare_input(input: str) -> str: 
    """
    Prepares the input for folding. 
    A value error is raised if any invalid letters are found in the sequence.  
    The prepared sequence is defined as a global variable

    args: 
        input: RNA sequence 
    """
    
    global sequence
    #Strips white space, changes to upper case and replaces T with U
    sequence = input.strip().upper().replace('T', 'U')

    allowed = set(['A', 'C', 'G', 'U', 'N'])
    if any(letter not in allowed for letter in sequence): 
        raise ValueError("Invalid letters found in the sequence")
    
    return sequence

def read_parameters(file_loop: str, file_stacking: str) -> tuple[pd.array, pd.array]: 
    """
    Read parameters from two .csv files:
    The .csv files are converted to Pandas tables and defined as the global variable 'parameters'

    Args: 
        file_loop: csv file containing the parameters for different loop types and lengths
        file_stacking: csv file containing the parameters for base paring between different nucleotides
    """
    try:
        loops = pd.read_csv(file_loop) 
        stacking = pd.read_csv(file_stacking, index_col=0)
    except FileNotFoundError:
        raise FileNotFoundError("One or both parameter files not found")
    except pd.errors.EmptyDataError:
        raise ValueError("One or both parameter files are empty or in an unexpected file format")
    
    global parameters
    parameters = (loops, stacking)
    
    return parameters

def write_dbn(name: str, sequence: str, fold: str, outfile: TextIOWrapper) -> None:
    """
    Writes sequence and fold to dbn file

    args: 
        name: name or ID of the sequence 
        sequence: RNA sequence
        fold: dot bracket structure
        outfile: an opened file to write to
    """
    outfile.write(f"#Name: {name}\n")
    outfile.write(f"#Length: {len(sequence)}\n")
    outfile.write(sequence + "\n")
    outfile.write(fold + "\n")


def parse_asymmetry_parameters(input_string: str) -> tuple[list, int]: 
    """
    Function for parsing the parameters for the penalty function for asymmetric loops that is given trough the command line
    Returns a list and an integer

    Args: 
        - input_string: string that represents a list and an integer. Example: "[0.4, 0.3, 0.2, 0.1]; 6]
    """
    try: 
        #Split into the f and penalty max parts
        arguments = input_string.split(';') 

        #Evaluate the input safely and evaluate literals. Can handle lists. 
        f = ast.literal_eval(arguments[0])
        assert type(f) == list

        penalty_max = int(arguments[1])

        return f, penalty_max
    
    except Exception as e: 
        raise argparse.ArgumentTypeError(f"Invalid input: {input_string}. Error: {str(e)}")

def declare_global_variable(b_stacking: bool = False, closing: bool = False, asymmetry: bool = False) -> None: 
    """
    Declares the global variables used troughout all the other functions. 

    Args: 
        - b_stacking: Is stacking retained for bulge loops of size 1 [True/False] (default = False)
        - closing: Is closing penalty added for GU/UG and AU/UA base pairs that closses interior loops [True/False] (default = False)
        - asymmetry: Is a penalty added for asymmetric interior loops [True/False] (default = False)
    """
    global basepairs, bulge_stacking, closing_penalty, asymmetry_penalty

    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}
    bulge_stacking = b_stacking
    closing_penalty = closing
    asymmetry_penalty = asymmetry



####### FOLD FUNCTIONS #######
def make_asymmetric_penalty(f: list, penalty_max: int) -> callable: 
    """
    f has to be a list. In articles writen as f(1), f(2) and so forth

    The asymmetry function is used to calculate a penalty to internal loops that are asymmetric. 
    This penalty does not exists in the orginial paper, but is added later

    This functions returns a function that uses the given parameters to calculate the penalty for asymmetric loops of given size

    args: 
        - f: list that gives the values from f(0) to f(m_max). 
        - penalty_max: that maximum penalty that can be addded to ad asymmetric interior loop
    """

    M_max = len(f)

    def asymmetry_func(i: int, ip: int, j: int, jp: int) -> float:
        """
        Calculates the penalty to add to the asymmetric interior loop enclosed by the base pairs i and j and i' and p' 
        """
        N1 = (ip-i-1)
        N2 =(j-jp-1)
        N = abs(N1-N2)
        M = min(M_max, N1, N2)-1
        penalty = min(penalty_max, N*f[M]) 
        return penalty
    
    global asymmetric_penalty_function
    asymmetric_penalty_function = asymmetry_func

    return asymmetric_penalty_function

def loop_greater_10(loop_type: str, length: int) -> float:
    """
    Calculates the energy parameters for loops with a size greater than 10 
    The parameter is calculated as described in 'Improved predictions of secondary structures for RNA'

    Args: 
        loop_type: type of loop to calculate energy for (HL, IL or BL)
        length: length of the loop
    """
    R = 0.001987 #In kcal/(mol*K)
    T = 310.15 #In K
    G_max = parameters[0].at[10, loop_type]

    G = G_max + 1.75*R*T*math.log(length/10)

    return G

### LOOP ENERGIES ###
def stacking(i: int, j: int, V: np.array) -> float: 
    """
    Find the energy parameter for basepairing of Sij and Si+1j-1, which results in basepair stacking
    If Si+1 and Sj+1 cannot basepair the energy is infinity
    Allows for Watson-Crick basepairs and wobble basepairs
    """

    prev_bp = sequence[i+1] + sequence[j-1]   
    
    #If previous bases can form a base pair stacking is possible
    if prev_bp in basepairs: 
        current_bp = sequence[i] + sequence[j]
        energy = round(parameters[1].at[current_bp, prev_bp] + V[i+1, j-1], 5)
    
    else: 
        energy = float('inf')
    
    return energy

def bulge_loop_3end(i: int, j: int, V: np.array) -> tuple[float, int]: 
    """
    Find the energy parameter of introducing a bulge loop on the 3' end of the strand 
    """
    energy = float('inf')
    ij = None

    #Try all sizes of bulge loop and save the one that gives the lowest energy
    for jp in range(i+2,j-1):  
        bp = sequence[i+1]+sequence[jp]
        if bp in basepairs: 
            size = j-jp-1
            if size <= 10:
               BL_energy = parameters[0].at[size, "BL"] + V[i+1, jp]
               if size == 1 and bulge_stacking: #Add stacking parameter if stacking for bulge loops is retained
                   BL_energy += parameters[1].at[(sequence[i]+sequence[j]), bp]
            else: 
                BL_energy = loop_greater_10("BL", size) + V[i+1, jp]
            
            if BL_energy < energy: 
                energy = BL_energy
                ij = jp
    
    return round(energy, 5), ij

def bulge_loop_5end(i: int, j: int, V: np.array) -> tuple[float, int]:
    """
    Find the energy parameter of introducing a bulge loop on the 5' end of the strand 
    """
    energy = float('inf')
    ij = None
    
    #Try all sizes of bulge loop and save the one that gives the lowest energy
    for ip in range(i+2,j-1):  
        bp = sequence[ip]+sequence[j-1]
        if bp in basepairs: 
            size = ip-i-1
            if size <= 10:
                BL_energy = parameters[0].at[size, "BL"] + V[ip, j-1]
                if size == 1 and bulge_stacking: #Add stacking parameter if stacking for bulge loops is retained
                    BL_energy += parameters[1].at[(sequence[i]+sequence[j]), bp] 
            else: 
                BL_energy = loop_greater_10("BL", size) + V[ip, j-1]

            if BL_energy < energy: 
                energy = BL_energy
                ij = ip

    return round(energy, 5), ij

def interior_loop(i: int, j: int, V: np.array) -> tuple[float, tuple[int, int]]: 
    """
    Find the energy parameter of adding a interior loop
    """
    
    energy = float('inf')
    ij = None

    for ip in range(i+2, j-2): #Try loop of any size between i and i'
        for jp in range(ip+3, j-1): #Try loop of any size between j and j'
            bp_prime = sequence[ip] + sequence[jp]
            if bp_prime in basepairs:
                size = (ip-i-1)+(j-jp-1)

                IL_energy = (parameters[0].at[size, "IL"] + V[ip, jp]) if size <= 10 else (loop_greater_10("IL", size) + V[ip, jp])
                
                #Add penalty to energy if loop is asymmetric
                if asymmetry_penalty and ((ip-i-1) != (j-jp-1)): 
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
                    energy = IL_energy
                    ij = (ip, jp)
    
    return round(energy, 5), ij

def find_E1(i: int, j: int) -> float:
    """
    E1 are the energy of base pairing between Si and Sj with one internal edge (hairpin loop) 
    """
    size = j-i-1    

    energy = parameters[0].at[size,"HL"] if size <= 10 else loop_greater_10("HL", size)

    return round(energy, 5)

def find_E2(i: int, j: int, V: np.array) -> float: 
    """
    E2 is the energy of basepairing between i and j and i' and j' resulting in two internal edges (stacking, bulge loop or internal loop)
    i<i'<j'<j
    Returns the minimum of the 3 options  
    """
    energy = min(stacking(i, j, V), 
                 bulge_loop_3end(i, j, V)[0], 
                 bulge_loop_5end(i, j, V)[0], 
                 interior_loop(i, j, V)[0])
    return energy

def find_E3(i: int, j: int, W: np.array) -> float: 
    """
    E3 is the energy of a structure that contains more than two internal edges (bifurcating loop)
    The energy is the energy of the sum of the substructures 
    i+1<i'<j-2
    """
    energy = float('inf')
    ij = None

    #Try all combinations of substructure and save the one that gives the lowest energy
    for ip in range(i+2, j-2):  
        loop_energy = W[i+1, ip] + W[ip+1, j-1]
        if loop_energy < energy: 
            energy = round(loop_energy, 5)
            ij = (ip, ip+1)
    return energy, ij

def find_E4(i: int, j: int, W: np.array) -> tuple[float, tuple[int, int]]: 
    """
    E4 is the energy when i and j are both in base pairs, but not with each other. 
    It find the minimum of combinations of two possible subsequences containing i and j
    """
    energy = float('inf')
    ij = None

    for ip in range(i+1, j-1): 
        subsequence_energy = W[i, ip] + W[ip+1, j]
        
        if subsequence_energy < energy: 
            energy = round(subsequence_energy, 5)
            ij = (ip, ip+1)

    return energy, ij

def penta_nucleotides(W: np.array, V: np.array) -> None:
    """
    Fills out the first entries in the matrices V and W 
    The shortest possible subsequences are of length 5 and can only form hairpin loops of size 3 if i and j basepair
    """
    N = len(sequence)

    for i in range(0, N-4): 
        j = i+4
        bp = sequence[i]+sequence[j]
        if bp not in basepairs:
            V[i,j] = W[i,j ]= float('inf')
        else: 
            V[i,j] = W[i,j] = parameters[0].at[3, "HL"] 

### FILL V AND W ###
def compute_V(i: int, j: int, W: np.array, V: np.array) -> None: 
    """
    Computes the minimization over E1, E2 and E3, which will give the value at V[i,j]
    """

    if sequence[i] + sequence[j] in basepairs:
        v = min(find_E1(i, j), 
                find_E2(i, j, V), 
                find_E3(i, j, W)[0])

    else: 
        v = float('inf')

    V[i, j] = v

def compute_W(i: int, j: int, W: np.array, V: np.array) -> None:
    """
    Computes the minimization over possibilities for W and fills out the entry at W[i,j]
     Possibilities are: 
    - i or j in a structure (W[i+1, j] or W[i, j-1])
    - i and j basepair with each other (V[i,j])
    - i and j both base pair but not with each other (E4)
    """
    w = min(W[i+1,j], W[i,j-1], V[i,j], find_E4(i, j, W)[0])

    W[i,j] = w


def fold_rna() -> tuple[np.array, np.array]: 
    """
    Fills out the W and V matrices to find the fold that gives the minimum free energy
    Follows Mfold as desribed by M. Zuker

    The V matrix contains the minimum free energy for the subsequences i and j, if i and j has to form a pair. 
    If i and j are not able to basepair the energy will be infinity (not a possible structure)

    The W matrix contains the minimum free energy for the subsequences i and j where base pairing between i and j is not nessecary.
    """
    N = len(sequence)
    W, V = np.full([N, N], float('inf')), np.full([N, N], float('inf'))


    #Fills out the table with all posible penta nucleotide subsequences
    # Penta nucleotides are the base cases. If subsequences are shorter they cannot be folded
    penta_nucleotides(W, V) 

    for l in range(5, N): #Computes the best score for all subsequences that are longer than 5 nucleotides with increasing length
        for i in range(0, N-5): 
            j = i+l
            if j < N: 
                compute_V(i, j, W, V) 
                compute_W(i, j, W, V)

    return W, V

def find_optimal(W: np.array) -> float: 
    """
    Find the final energy of the folded RNA
    """
    return W[0, -1]

### BACTRACKING ### 
def trace_V(i: int, j: int, W: np.array, V: np.array, dotbracket: list) -> None: 
    """
    Traces backwards trough the V matrix recursively to find the secondary structure
    """
    if V[i,j] == find_E1(i, j): 
        dotbracket[i], dotbracket[j] = '(', ')'
        for n in range(i+1, j): 
            dotbracket[n] = '.'
    
    elif V[i,j] == stacking(i, j, V): 
        dotbracket[i], dotbracket[j] = '(', ')'
        trace_V(i+1, j-1, W, V, dotbracket)
    
    elif V[i,j] == bulge_loop_3end(i, j, V)[0]: 
        jp = bulge_loop_3end(i, j, V)[1]
        dotbracket[i], dotbracket[j] = '(', ')'
        for n in range(jp, j): 
            dotbracket[n] = '.'
        trace_V(i+1, jp, W, V)
    
    elif V[i,j] == bulge_loop_5end(i, j, V)[0]: 
        ip = bulge_loop_5end(i, j, V)[1]
        dotbracket[i], dotbracket[j] = '(', ')'
        for n in range(i+1, ip): 
            dotbracket[n] = '.'
        trace_V(ip, j-1, W, V, dotbracket)
    
    elif V[i,j] == interior_loop(i, j, V)[0]:
        ij = interior_loop(i, j, V)[1]
        dotbracket[i], dotbracket[j] = '(', ')' 
        for n in range(i+1, ij[0]): 
            dotbracket[n] = '.'
        for n in range(ij[1]+1, j): 
            dotbracket[n] = '.'
        trace_V(ij[0], ij[1], W, V, dotbracket)
    
    elif V[i, j] == find_E3(i, j, W)[0]: 
        ij = find_E3(i, j, W)[1]
        dotbracket[i], dotbracket[j] = '(', ')' 
        trace_W(i+1, ij[0], W, V, dotbracket), trace_W(ij[1], j-1, W, V, dotbracket)

def trace_W(i: int, j: int, W: np.array, V: np.array, dotbracket: list) -> None: 
    """
    Traces backwards trough the W matrix recursively to find the secondary structure
    """
    if W[i,j] == W[i+1, j]: 
        dotbracket[i] = '.'
        trace_W(i+1, j, W, V, dotbracket)

    elif W[i,j] == W[i, j-1]: 
        dotbracket[j] = '.'
        trace_W(i, j-1, W, V, dotbracket)

    elif W[i, j] == V[i, j]: 
        trace_V(i, j, W, V, dotbracket)

    elif W[i,j] == find_E4(i, j, W)[0]: 
        ij = find_E4(i,j,W)[1] 
        trace_W(i, ij[0], W, V, dotbracket), trace_W(ij[1], j, W, V, dotbracket)



def backtrack(W: np.array, V: np.array) -> str: 
    """
    Backtracks trough the W, V matrices to find the final fold
    Returns the fold as a dotbracket structure
    """
    #Allocate the dot bracket structure
    dotbracket =  ['?' for x in range(W.shape[0])]
    
    j = W.shape[0]-1
    i = 0
    
    #Fill out db
    trace_W(i, j, W, V,dotbracket)

    return "".join(dotbracket)



def main() -> None: 
    """
    Running the program to predict the secondary structure of an RNA sequence
    
    args:
        - i: input given as sequence in command line 
        - f: input given as fasta file containing one seqeuence
        - lp: specifies which set of parameters is used [1988, 1989] (default = 1989)
        - b: used to retain base pair stacking for bulge loops of size one
        - a: used to add penalty for asymmetric interior loops
        - c: used to add penalty for GU/UG and AU/UA base pairs that closes interior loops
        - A: used to specify the parameters to use to calculate asymmetric penalty ["f; penalty_max"] (default = "[0.4, 0.3, 0.2, 0.1], 3")
        - lf and sf: used to add costum parameter files for loops and stacking respectively
        - o: name of output file for dot bracket structure without extension (default = stdout)
    """
    #Setting up the option parsing using the argparse module
    argparser = argparse.ArgumentParser(prog='Mfold',
        description="""Version of Mfold implemented by Maria Eskerod.
                     Implemented as a part of a project in bioinformatics, fall 2023.
                     Bioinformatic Research Center at Aarhus University""", 
        epilog="""Penalty for asymmetry for interior loops: A penalty is calculated as min(max penalty, N*f(M)) 
        where M = min(max_M, N1, N2) with N1 and N2 being the number of unpaired bases on each strand and N = abs(N1 - N2)""")
    #Adding arguments
    #Input can either be provided in a file or in stdin
    argparser.add_argument('-i', '--input', metavar='', help = 'input provided as sequence in command line') 
    argparser.add_argument('-f', '--file', type=argparse.FileType('r'), metavar='', help = 'input provided as fasta file')
    #Used to modify which version of Mfold is used
    argparser.add_argument('-lp', '--loop_parameters', type=str, choices=['1988', '1989'], default='1989', help='used to specify what set of loop parameters should be used. The options are 1988 and 1989 (default = 1988)')
    argparser.add_argument('-b', '--bulge_stacking', action='store_true',  help='if -b/--bulge_stacking is used stacking is retained for bulge loops of size 1 and the energy parameter for stacking will be added to the loop energy')
    argparser.add_argument('-a', '--asymmetric', action='store_true', help='if -a/--asymmetric is used a penalty will be added for asymmetric interior loops (see below)')
    argparser.add_argument('-c', '--closing_penalty', action='store_true', help='if -c/--closing_penalty is used a penalty will be added for GU/UG and AU/UA base pairs that closes interior loops')
    argparser.add_argument('-A', '--asymmetry_parameters', type=parse_asymmetry_parameters, default=[[0.4, 0.3, 0.2, 0.1], 3], metavar='', help='parameters used to calculate the penalty to add for asymmetric interior loops. Must be supplied as "f; penalty_max" (default = "[0.4, 0.3, 0.2, 0.1], 3") (see further description below). Only used if -a is also used')
    argparser.add_argument('-lf', '--loop_file', metavar= '', help='user defined parameter file for loops. Must contain energy for hairpin (HL), interior (IL) and bulge loops (BL) for sizes up to 10 in a csv file.')
    argparser.add_argument('-sf', '--stacking_file', metavar='', help='User defined parameter file for base pair stacking. Must contain energy for all possible base pairs as a matrix in a csv file.')
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('-o', '--outfile', metavar='', default=sys.stdout, help='name for output file without extension')

    args = argparser.parse_args()

    #Read sequence
    if args.input: 
        prepare_input(args.input)
        name = "user inputted sequence"

    if args.file: 
        raw_sequence, name = read_fasta(args.file)
        prepare_input(raw_sequence)
    
    if not sequence: 
        raise ValueError("No valid input sequence provided.")
    
    #Declare the global parameters that defines the version of Mfold
    declare_global_variable(args.bulge_stacking, args.closing_penalty, args.asymmetric)

    #Make penalty function for asymmetric loops
    f, penalty_max = args.asymmetry_parameters
    make_asymmetric_penalty(f, penalty_max)


    #Find absolute path of script, to be able to open csv files whereever the script is called from
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    #Define the parameters to use for folding
    if args.loop_file: 
        loop_file = args.loop_file
    elif args.loop_parameters == '1988':
        loop_file = os.path.join(script_dir, 'parameters', 'loop_1988.csv')
    elif args.loop_parameters == '1989':
        loop_file = os.path.join(script_dir, 'parameters', 'loop_1989.csv')


    if args.stacking_file:
        stacking_file = args.stacking_file
    else: 
        stacking_file = os.path.join(script_dir, 'parameters', 'stacking_1988.csv')

    read_parameters(loop_file, stacking_file)

    
    
    print(f"Fold {name}\n")
    start_time = time.time()

    W, V = fold_rna()
    energy = find_optimal(W)
    fold = backtrack(W, V) 

    #Write to outfile
    if args.outfile == sys.stdout: 
        args.outfile.write(fold + "\n\n")
    
    else: 
        outfile = args.outfile + ".dbn"
        with open(outfile, 'w') as file:
            write_dbn(name, sequence, fold, file)

    print("Finished in {} second".format(round(time.time() - start_time, 2)))
    print(f"Energy of optimal fold is {energy} kcal/mol\n") 

if __name__ == '__main__': 
    main()