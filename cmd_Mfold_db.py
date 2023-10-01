from io import TextIOWrapper
import argparse, subprocess, os

def run_Mfold(filename): 
    """
    """
    command = f'mfold SEQ={filename}'

    result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    result.wait()

def read_ct(file) -> tuple():
    """
    Takes a .ct file and returns the sequence as a string and a list of base pairs
    """
    length = 0
    pairs = []

    with open(file, 'r') as f:
        lines = [line.split() for line in f.readlines()]

    #Remove header - if any
    header_lines = 0
    for line in lines: 
        if line[0] == '1': 
                break
        else: 
            header_lines += 1

    lines = lines[header_lines:]
    
    for line in lines: 
        length += 1
        if line[4] != '0': 
            pairs.append((int(line[0])-1, int(line[4])-1)) #The files start indexing from 1

    return length, pairs

def remove_files(filename): 
    """
    """
    target, _ = os.path.splitext(os.path.basename(filename))
    cdir = os.getcwd()
    for name in os.listdir(cdir):
         if os.path.isfile(name):
              if name.startswith(target):
                   os.remove(os.path.join(cdir, name))

def ct_to_db(length, pairs):
    db = ['.' for n in range(length)] 

    for pair in pairs: 
        if pair[0] < pair[1]: 
             db[pair[0]] = '('
             db[pair[1]] = ')' 
     
    return "".join(db)

def db_to_file(db, filename, name): 
    with open(filename, 'w') as f: 
        f.write(f">{name}\n")
        f.write(f"{db}\n")

def main(): 
    argparser = argparse.ArgumentParser(description="" )
    #Adding arguments
    #TODO - Add description/help for command line options
    #Input can either be provided in a file or in stdin
    argparser.add_argument('filename') 
    #Setting up output. Writes to specified outfile or stdout
    argparser.add_argument('outputdir')

    args = argparser.parse_args()

    file = args.filename
    run_Mfold(file)

    basename = os.path.splitext(os.path.basename(file))[0]
    ct_file = basename + '.ct'

    length, pairs = read_ct(ct_file) 

    remove_files(file)

    db = ct_to_db(length, pairs)

    outfile = os.path.join(args.outputdir, basename + '.txt')

    db_to_file(db, outfile, basename)

if __name__ == '__main__': 
     main()