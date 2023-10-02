import argparse, time

from main import(read_fasta, prepare_input)

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
        name = "User inputted sequence"

    if args.file: 
        sequence, name = read_fasta(args.file)
        sequence = prepare_input(sequence)
    
    if not sequence: 
        raise ValueError("No valid input sequence provided.")

    print(f"Fold {name}")
    start_time = time.time()

    #FIXME - Add version of Mfold which depends on an argument

    print("Finished in {} second".format(round(time.time() - start_time, 2)))
    print(f"Energy of optimal fold is {"energy"} kcal/mol\n") #NOTE - Change energy

    #Write to outfile
    args.outfile.write(f">{name}\n")
    args.outfile.write(f"{"fold"}\n") #NOTE - Change fold


if __name__ == '__main__': 
    main()