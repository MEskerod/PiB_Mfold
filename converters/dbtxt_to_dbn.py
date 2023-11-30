import os
from Bio import SeqIO


def get_path_list(dir_name: str) -> list[str]: 
    """
    Makes a list of all the files in the given directory. 
    Returns a list of "dir_name/file_name"
    """
    file_names = [os.path.join(dir_name, name) for name in os.listdir(dir_name)]
    return file_names

def read_fasta(input: str) -> tuple[str, str]:
    """
    Reads in a FASTA-file and returns the sequence
    If there is more than one sequence in the FASTA file it gives an error
    """
    with open(input, 'r') as f:
        records = list(SeqIO.parse(f, 'fasta'))
    
    if len(records) > 1: 
        raise ValueError("FASTA file contains more than one sequence")

    return str(records[0].seq), records[0].id

def read_dbtxt(file_name: str) -> str: 
    """
    Read a dot bracket structure from a txt file
    Header lines starts with >
    """
    lines = []
    with open(file_name, 'r') as f: 
        for line in f:
            if line.startswith('>'): 
                continue
            lines.append(line.strip())
    return "".join(lines) 


def write_dbn(name: str, sequence: str, fold: str, outfile: str) -> None:
    """
    Writes name, sequence and structure to a dbn file.
    Header lines are denoted by #
    """
    with open(outfile, 'w') as f:
        f.write(f"#Name: {name}\n")
        f.write(f"#Length: {len(sequence)}\n")
        f.write(sequence + "\n")
        f.write(fold + "\n")

def test_matching_filelist(filelist1: list[str], filelist2: list[str]) -> None:
    """
    Checks that the the two file lists contains the same files
    """
    assert len(filelist1) == len(filelist2)

    for n in range(len(filelist1)): 
        name1 = os.path.splitext(os.path.basename(filelist1[n]))[0]
        name2 = os.path.splitext(os.path.basename(filelist2[n]))[0]
        assert name1 == name2


def main() -> None: 
    """
    Takes the sequences and structures, merges them and saves them as a dbn file
    """
    seuence_dir = "../sequences"
    dbtxt_dir = "../true_structures"

    ouput_dir = "../all_structures"

    os.makedirs(ouput_dir, exist_ok=True)

    sequence_filelist = get_path_list(seuence_dir)
    dbtxt_filelist = get_path_list(dbtxt_dir)

    test_matching_filelist(sequence_filelist, dbtxt_filelist)

    for n in range(len(sequence_filelist)):   
        sequence, id = read_fasta(sequence_filelist[n])
        db = read_dbtxt(dbtxt_filelist[n])
        basename = os.path.splitext(os.path.basename(sequence_filelist[n]))[0]
        output_file = os.path.join(ouput_dir, basename + '.dbn')
        assert len(sequence) == len(db)
        write_dbn(id, sequence, db, output_file)

if __name__ == '__main__': 
    main()
