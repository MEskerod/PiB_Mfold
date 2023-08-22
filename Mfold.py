import argparse, sys, fileinput
import numpy as np

def prepare_input(input: list[str]) -> str: 
    """
    Removes name of sequence if input is fasta 
    Strips any whitespace
    Makes sure that all letters are in upper case 
    """
    input = "".join([line.strip().upper() for line in input if line[0].isalpha()])
    return input

def read_parameters(file): 
    return

def fold_rna(): 
    return

def main() -> None: 
    """
    Running the program
    Takes following arguments: 
    input - Can be a file. If no file is given it's possible to write the sequence in stdin after running the program and end by pressing ctrl+d
    ouput - Can be a file specified by using the flag -o or a the default stdout
    """
    #Setting up the option parsing using the argparse module
    argparser = argparse.ArgumentParser(
        description="" )
    #Adding arguments
    #TODO - Set up arguments
    #Input can either be provided in a file or in stdin
    argparser.add_argument('files', nargs='*')
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('-o', '--outfile', metavar='output', type=argparse.FileType('w'), default=sys.stdout)

    args = argparser.parse_args()

    sequence = []
    for line in fileinput.input(files=args.files): 
        sequence.append(line)
    
    sequence = prepare_input(sequence)

    parameters = read_parameters("parameters.txt")
    
    print(sequence)

    

    return

if __name__ == '__main__': 
    main()