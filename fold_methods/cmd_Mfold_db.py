from io import TextIOWrapper
from Bio import SeqIO
import argparse, subprocess, os

def run_Mfold(filename: str) -> None: 
    """
    Runs the version of Mfold downloaded from unafold website
    """
    command = f'mfold SEQ={filename}'

    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    result.wait()

def read_ct(file:str) -> tuple:
    """  
    Takes a .ct file and returns the sequence as a string and a list of base pairs
    """
    pairs = []
    sequence = ""

    with open(file, 'r') as f:
        lines = [line.split() for line in f.readlines()]

    length = int(lines[0][0])
    
    for n in range(1, length+1):
        sequence += lines[n][1]
        if lines[n][4] != '0': 
            pairs.append((int(lines[n][0])-1, int(lines[n][4])-1)) #The files start indexing from 1

    return sequence, pairs

def sequence_pair_from_fasta(file: str) -> tuple: 
    """
    Reads in a FASTA-file and returns the sequence
    If there is more than one sequence in the FASTA file it gives an error
    It returns the sequence and an empty list representing the possible pairs
    """
    records = list(SeqIO.parse(file, 'fasta'))
    
    if len(records) > 1: 
        raise ValueError("FASTA file contains more than one sequence")
    
    return str(records[0].seq), []

def remove_files(filename: str) -> None: 
    """
    Remove file filename
    """
    target, _ = os.path.splitext(os.path.basename(filename))
    cdir = os.getcwd()
    for name in os.listdir(cdir):
         if os.path.isfile(name):
              if name.startswith(target):
                   os.remove(os.path.join(cdir, name))

def ct_to_db(length: int, pairs: list[tuple]):
    """
    Takes a list of base pairs and converts it to a db file
    Cannot handle pseudo knots
    """
    db = ['.' for n in range(length)] 

    for pair in pairs: 
        if pair[0] < pair[1]: 
             db[pair[0]] = '('
             db[pair[1]] = ')' 
     
    return "".join(db)

def db_to_file(sequence: str, db: str, filename: str, name: str) -> None: 
    """
    Writes a sequence and structure to a dbn file 
    Header lines a denoted by #
    """
    with open(filename, 'w') as f: 
        f.write(f"#Name: {name}\n")
        f.write(f"#Length: {len(sequence)}\n")
        f.write(sequence + "\n")
        f.write(db + "\n")


def main(): 
    argparser = argparse.ArgumentParser(description="Program to run the version of Mfold from unafold website (anno fall 2023) and save the structure as dbn file " )
    #Adding arguments
    #Input can either be provided in a file or in stdin
    argparser.add_argument('filename', help="Name of input file") 
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('outputdir', help = "Name of directory where file is outputted")

    args = argparser.parse_args()

    file = args.filename
    run_Mfold(file)

    basename = os.path.splitext(os.path.basename(file))[0]
    ct_file = basename + '.ct'

    #If Mfold isn't able to generate a structure the string is rendered unfolded
    if os.path.exists(ct_file):
        sequence, pairs = read_ct(ct_file) 
    else: 
        sequence, pairs = sequence_pair_from_fasta(file)

    remove_files(file)

    db = ct_to_db(len(sequence), pairs)

    outfile = os.path.join(args.outputdir, basename + '.dbn')

    db_to_file(sequence, db, outfile, basename)

if __name__ == '__main__': 
     main()