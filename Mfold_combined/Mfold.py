import ast, argparse, sys, time, math
import numpy as np
import pandas as pd
from Bio import SeqIO
from io import TextIOWrapper

####### HELP FUNCTIONS #######
def read_fasta(input) -> str:
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
    
    global sequence
    sequence = input.strip().upper().replace('T', 'U')

    allowed = set(['A', 'C', 'G', 'U', 'N'])
    if any(letter not in allowed for letter in sequence): 
        raise ValueError("Invalid letters found in the sequence")
    
    return sequence

def read_parameters(file_loop, file_stacking): 
    """
    Read parameters from to .csv files:
    - One containing the loop parameters 
    - A second one containg the parameters for different types og base pairing 
    The .csv files are converted to Pandas tables and returned
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

def write_dbn(name, sequence, fold, outfile: TextIOWrapper):
    """
    """
    outfile.write(f"#Name: {name}\n")
    outfile.write(f"#Length: {len(sequence)}\n")
    outfile.write(sequence + "\n")
    outfile.write(fold + "\n")


def parse_asymmetry_parameters(input_string: str): 
    """
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

def declare_global_variable(b_stacking = False, closing = False, asymmetry = False): 
    global basepairs, bulge_stacking, closing_penalty, asymmetry_penalty

    basepairs = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}
    bulge_stacking = b_stacking
    closing_penalty = closing
    asymmetry_penalty = asymmetry



####### FOLD FUNCTIONS #######
def make_asymmetric_penalty(f, penalty_max): 
    """
    f has to be a list. In articles writen as f(1), f(2) and so forth

    The asymmetry function is used to calculate a penalty to internal loops that are asymmetric. 
    This penalty does not exists in the orginial paper, but is added later

    This functions returns a function that uses the given parameters to calculate the penalty for asymmetric loops of given size
    """

    M_max = len(f)

    def asymmetry_func(i, ip, j, jp):
        N1 = (ip-i-1)
        N2 =(j-jp-1)
        N = abs(N1-N2)
        M = min(M_max, N1, N2)-1
        penalty = min(penalty_max, N*f[M]) 
        return penalty
    
    global asymmetric_penalty_function
    asymmetric_penalty_function = asymmetry_func

    return asymmetric_penalty_function

def loop_greater_10(loop_type, length):
    """
    Calculates the energy parameters for loops with a size greater than 10 
    The parameter is calculated as described in 'Improved predictions of secondary structures for RNA'

    The function returns the energi parameter for a loop of a given type and size. 
    """
    R = 0.001987 #In kcal/(mol*K)
    T = 310.15 #In K
    G_max = parameters[0].at[10, loop_type]

    G = G_max + 1.75*R*T*math.log(length/10)

    return G

### LOOP ENERGIES ###
def stacking(i, j, V): 
    """
    Find the energy parameter for basepairing of Sij and Si+1j-1, which results in basepair stacking
    If Si+1 and Sj+1 cannot basepair the energy is infinity
    Allows for Watson-Crick basepairs and wobble basepairs
    """

    prev_bp = sequence[i+1] + sequence[j-1]   
    
    if prev_bp in basepairs: 
        current_bp = sequence[i] + sequence[j]
        energy = round(parameters[1].at[current_bp, prev_bp] + V[i+1, j-1], 5)
    else: 
        energy = float('inf')
    return energy

def bulge_loop_3end(i, j, V): 
    """
    Find the energy parameter of introducing a bulge loop on the 3' end. 
    If the size of the bulge loop is 1 and stacking=True, stacking on each side of the loop is retained and stacking parameter is added
    Is able to handle loops of any size
    """
    energy = float('inf')
    ij = None

    #Bulge on 3' end
    for jp in range(i+2,j-1):  
        bp = sequence[i+1]+sequence[jp]
        if bp in basepairs: 
            size = j-jp-1
            if size <= 10:
               BL_energy = parameters[0].at[size, "BL"] + V[i+1, jp]
               if size == 1 and bulge_stacking: 
                   BL_energy += parameters[1].at[(sequence[i]+sequence[j]), bp]
            else: 
                BL_energy = loop_greater_10("BL", size) + V[i+1, jp]
            
            if BL_energy < energy: 
                energy = round(BL_energy, 5)
                ij = jp
    
    return energy, ij

def bulge_loop_5end(i, j, V):
    """
    Find the energy parameter of introducing a bulge loop on the 5' end. 
    If the size of the bulge loop is 1 and stacking=True, stacking on each side of the loop is retained and stacking parameter is added
    Is able to handle loops of any size
    """
    energy = float('inf')
    ij = None
    
    #Bulge on 5' end
    for ip in range(i+2,j-1):  
        bp = sequence[ip]+sequence[j-1]
        if bp in basepairs: 
            size = ip-i-1
            if size <= 10:
                BL_energy = parameters[0].at[size, "BL"] + V[ip, j-1]
                if size == 1 and bulge_stacking: 
                    BL_energy += parameters[1].at[(sequence[i]+sequence[j]), bp] 
            else: 
                BL_energy = loop_greater_10("BL", size) + V[ip, j-1]

            if BL_energy < energy: 
                energy = round(BL_energy, 5)
                ij = ip

    return energy, ij

def interior_loop(i, j, V): 
    """
    Find the energy parameter of adding a interior loop. 
    Is able to handle loops of any size
    A penalty is added for asymmetric loops. If the penalty should not be added interior_loop should be called with asymmetry = False
    A penalty is added for interior loops closed by AU and GU base pairs
    """
    
    energy = float('inf')
    ij = None

    for ip in range(i+2, j-2): 
        for jp in range(ip+3, j-1): 
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
                    energy = round(IL_energy, 5)
                    ij = (ip, jp)
    
    return energy, ij

def find_E1(i, j):
    """
    E1 is when Si and Sj basepair and gives one internal edge 
    This gives a hairpin loop 
    The function is able to handle loops of any size
    """
    size = j-i-1    

    energy = parameters[0].at[size,"HL"] if size <= 10 else loop_greater_10("HL", size)

    return round(energy, 5)

def find_E2(i, j, V): 
    """
    E2 is when Si and Sj contains two internal edged 
    Contains an edge between Si and Sj and an edge between Si' and Sj' 
    i<i'<j'<j
    Can be stacking, bulge loops or internal loops
    Returns the minimum of the 3 options. 
    """
    energy = min(stacking(i, j, V), 
                 bulge_loop_3end(i, j, V)[0], 
                 bulge_loop_5end(i, j, V)[0], 
                 interior_loop(i, j, V)[0])
    return energy

def find_E3(i, j, W): 
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
def penta_nucleotides(W, V):
    """
    Initiates the V and W matrices. 
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

def compute_V(i, j, W, V): 
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


def fold_rna(): 
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
    penta_nucleotides(W, V) 

    for l in range(5, N): #Computes the best score for all subsequences that are longer than 5 nucleotides
        for i in range(0, N-5): 
            j = i+l
            if j < N: 
                compute_V(i, j, W, V) 
                compute_W(i, j, W, V)

    return W, V

def find_optimal(W) -> float: 
    """
    Find the final energy of the folded RNA
    """
    return W[0, -1]

### BACTRACKING ### 
def trace_V(i, j, W, V, dotbracket): 
    """
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

def trace_W(i, j, W, V, dotbracket): 
    """
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



def backtrack(W, V): 
    """
    Backtracks trough the W, V matrices to find the final fold
    """
    dotbracket =  ['?' for x in range(W.shape[0])]
    
    j = W.shape[0]-1
    i = 0
    
    trace_W(i, j, W, V,dotbracket)

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
    #TODO - Add description/help for command line options
    #Input can either be provided in a file or in stdin
    argparser.add_argument('-i', '--input') 
    argparser.add_argument('-f', '--file', type=argparse.FileType('r'))
    argparser.add_argument('-lp', '--loop_parameters', type=str, choices=['1988', '1989'], default='1989')
    argparser.add_argument('-b', '--bulge_stacking', action='store_true')
    argparser.add_argument('-a', '--asymmetric', action='store_true')
    argparser.add_argument('-c', '--closing_penalty', action='store_true')
    argparser.add_argument('-A', '--asymmetry_parameters', type=parse_asymmetry_parameters, default=[[0.4, 0.3, 0.2, 0.1], 3])
    argparser.add_argument('-lf', '--loop_file')
    argparser.add_argument('-sf', '--stacking_file')
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('-o', '--outfile', metavar='output', default=sys.stdout)

    args = argparser.parse_args()

    if args.input: 
        prepare_input(args.input)
        name = "user inputted sequence"

    if args.file: 
        raw_sequence, name = read_fasta(args.file)
        prepare_input(raw_sequence)
    
    if not sequence: 
        raise ValueError("No valid input sequence provided.")
    
    declare_global_variable(args.bulge_stacking, args.closing_penalty, args.asymmetric)

    f, penalty_max = args.asymmetry_parameters

    #Change to find the right path whereever they are called from
    if args.loop_file: 
        loop_file = args.loop_file
    elif args.loop_parameters == '1988':
        loop_file = "../Mfold/parameters/loop_1988.csv"
    elif args.loop_parameters == '1989':
        loop_file = "../Mfold/parameters/loop_1989.csv"


    if args.stacking_file:
        stacking_file = args.stacking_file
    else: 
        stacking_file = "../Mfold/parameters/stacking_1988.csv"

    read_parameters(loop_file, stacking_file)

    make_asymmetric_penalty(f, penalty_max)

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