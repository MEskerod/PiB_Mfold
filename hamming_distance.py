"""
Script for calculating the structural Hamming distance between two dot bracket RNA secondary structures
The dot bracket structures should be in a .txt file 
Besides the dot bracket structure the only allowed in the file is a line with the name af the molecule: 
> "Name"
"""
import argparse

def read_dot_bracket(file):
    """
    Read the dot bracket structure in the text file into a string
    """
    dotbracket = "".join([line.strip() for line in file if not line.startswith(">")])
    return dotbracket


def calculate_distance(structure1, structure2): 
    """
    Calculates the Hamming distance between two sequences. 
    The distance is the sum of all the positions that does not have the same value
    """
    if len(structure1) != len(structure2):
        raise ValueError("Structures must have same length")
    
    distance = sum(s1 != s2 for s1, s2 in zip(structure1, structure2))

    return distance

def structural_hamming_distance(file1, file2): 
    """
    Calculated the structural Hamming Distance from two files
    """
    structure1 = read_dot_bracket(file1)
    structure2 = read_dot_bracket(file2)
    distance = calculate_distance(structure1, structure2)
    return distance


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
    argparser.add_argument('structure1', type=argparse.FileType('r'))
    argparser.add_argument('structure2', type=argparse.FileType('r'))
    #Setting up output. Writes to specified outfile or stdout

    args = argparser.parse_args()

    distance = structural_hamming_distance(args.structure1, args.structure2)

    print(distance)

    
if __name__ == '__main__': 
    main()


#print(read_dot_bracket("structures/67_URS00005D4795.txt"))