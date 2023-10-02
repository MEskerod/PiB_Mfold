import argparse, sys, time

from help_functions import(read_fasta, prepare_input, read_general_parameters, write_dbn)
from fold_functions import(loop_greater_10, make_asymmetric_penalty, fold_rna, find_optimal, backtrack) #TODO - Add backtrack when ready!

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
    #argparser.add_argument() #TODO - Add arguments for folding options and parameters
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('-o', '--outfile', metavar='output', default=sys.stdout)

    args = argparser.parse_args()

    if args.input: 
        sequence = prepare_input(args.input)
        name = "User inputted sequence"

    if args.file: 
        sequence, name = read_fasta(args.file)
        sequence = prepare_input(sequence)
    
    if not sequence: 
        raise ValueError("No valid input sequence provided.")
    
    #TODO - Change the below to match with the arguments passed by user
    bulge_stacking = True
    closing_penalty = True
    asymmetry_penalty = True

    #TODO - Change the below to match with the arguments passed by user
    loop_file = "parameters/loop_1989.csv"
    stacking_file = "parameters/stacking_1988.csv"

    parameters = read_general_parameters(loop_file, stacking_file)

    #FIXME - Add version of Mfold which depends on an argument
    #TODO - Change the below to match with the arguments passed by user
    asymmetric_penalty_function = make_asymmetric_penalty([0.4, 0.3, 0.2, 0.1], 3)

    print(f"Fold {name}\n")
    start_time = time.time()

    W, V = fold_rna(sequence, parameters, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)
    energy = find_optimal(W)
    fold = backtrack(W, V, parameters, sequence) #NOTE - May need to be updated if dangling ends are added

    #Write to outfile
    if args.outfile == sys.stdout: 
        args.outfile.write(fold + "\n\n")
    
    else: 
        outfile = args.outfile + ".dbn"
        with open(outfile, 'w') as f:
            write_dbn(name, sequence, fold, f)

    print("Finished in {} second".format(round(time.time() - start_time, 2)))
    print(f"Energy of optimal fold is {energy} kcal/mol\n") #NOTE - Change energy

if __name__ == '__main__': 
    main()