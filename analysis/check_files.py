import os
from Bio import SeqIO

from general import read_dbn_file, read_fasta


def get_path_list(dir_name): 
    """
    """
    file_names = [os.path.splitext(name)[0] for name in os.listdir(dir_name)]
    return file_names

#####################

def check_file_lists(file_list1, file_list2): 
    assert len(file_list1) == len(file_list2)

    for n in range(len(file_list2)): 
        try: 
            assert file_list1[n] == file_list2[n]
        except: 
            raise ValueError(f"ERROR: {file_list1[n]} is not equal to {file_list2[n]}")

def check_content(file_names, sequences_path, structures_path): 
    for n in range(len(file_names)): 
        file_name = file_names[n]
        length = int(file_name.split("_")[0])
        sequence, _ = read_fasta(os.path.join(sequences_path, file_name + ".fasta"))
        structure, _ = read_dbn_file(os.path.join(structures_path, file_name + ".dbn"))

        if (length != len(sequence)) & (length != len(structure)):
            raise ValueError(f"ERROR: {file_name} named wrong!")
        elif len(sequence) != len(structure): 
            raise ValueError(f"ERROR: Structure and sequence for {file_name} have different lengths")

def check_files(sequence_path, structures_path):
    filelist1 = get_path_list(sequence_path)
    filelist2 = get_path_list(structures_path)

    check_file_lists(filelist1, filelist2)

    check_content(filelist1, sequence_path, structures_path)
    return


def main(): 
    sequence_path = "../sequences"
    structures_path = "../structures/true_structures"

    check_files(sequence_path, structures_path)


if __name__ == "__main__": 
    main()