import argparse, sys, time

from help_functions import(read_fasta, prepare_input, read_general_parameters, write_dbn, parse_asymmetry_parameters)
from fold_functions import(loop_greater_10, make_asymmetric_penalty, fold_rna, find_optimal, backtrack) 

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
        sequence = prepare_input(args.input)
        name = "user inputted sequence"

    if args.file: 
        sequence, name = read_fasta(args.file)
        sequence = prepare_input(sequence)
    
    if not sequence: 
        raise ValueError("No valid input sequence provided.")
    
    bulge_stacking = args.bulge_stacking
    closing_penalty = args.closing_penalty
    asymmetry_penalty = args.asymmetric

    f, penalty_max = args.asymmetry_parameters

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

    parameters = read_general_parameters(loop_file, stacking_file)

    asymmetric_penalty_function = make_asymmetric_penalty(f, penalty_max)

    print(f"Fold {name}\n")
    start_time = time.time()

    W, V = fold_rna(sequence, parameters, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty)
    energy = find_optimal(W)
    fold = backtrack(W, V, parameters, sequence, asymmetric_penalty_function, bulge_stacking, closing_penalty, asymmetry_penalty) 

    #Write to outfile
    if args.outfile == sys.stdout: 
        args.outfile.write(fold + "\n\n")
    
    else: 
        outfile = args.outfile + ".dbn"
        with open(outfile, 'w') as f:
            write_dbn(name, sequence, fold, f)

    print("Finished in {} second".format(round(time.time() - start_time, 2)))
    print(f"Energy of optimal fold is {energy} kcal/mol\n") 

if __name__ == '__main__': 
    main()