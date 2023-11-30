import os

from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord

def read_dbn(file_name: str) -> tuple[str, str]:
    """
    Reads the dbn file and returns the sequence and id
    """
    lines = []
    with open(file_name, 'r') as f: 
        lines = f.readlines()

    sequence = lines[-2].strip()

    id = lines[0].split()[1]

    return sequence, id

def write_fasta(sequence: str, seq_id: str, file_name: str):
    """
    Writes a sequence to a fasta file with its id
    """
    seq = Seq(sequence)
    record =  SeqRecord(seq, id = seq_id)

    with open(file_name, "w") as f: 
        SeqIO.write(record, f, "fasta")

def make_file_name(sequence: str, input_file: str, output_dir: str):
    """
    Creates the file name of the new file as length_name.fasta, with the basename being the same as the original file
    """
    basename = os.path.splitext(os.path.basename(input_file))[0].replace('_', '')
    length = len(sequence)
    
    name = os.path.join(output_dir, "".join([str(length), '_', basename, '.fasta']))
    return name



def get_path_list(dir_name: str) -> list[str]: 
    """
    Makes a list of all the files in the given directory. 
    Returns a list of "dir_name/file_name"
    """
    file_names = [os.path.join(dir_name, name) for name in os.listdir(dir_name)]
    return file_names

def main() -> None: 
    path_list  = get_path_list("../dbn")
    output_dir = "../all_sequences"

    os.makedirs(output_dir, exist_ok=True)

    for file in path_list:
        sequence, id = read_dbn(file)
        file_name = make_file_name(sequence, file, output_dir)
        write_fasta(sequence, id, file_name)
    return


if __name__ == '__main__': 
    main()