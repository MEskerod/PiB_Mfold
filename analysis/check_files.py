import os, sys
from Bio import SeqIO

from general import read_dbn_file, read_fasta


def costum_sort(name: str) -> tuple[int, str, str]: 
    """
    Function to be used for sorting. 
    Ensures that the elements are sorted after the first element, the second element, and third element in the file name and not just the characters
    """
    components = name.split('_')
    return (int(components[0]), components[1], components[2])

def get_path_list(dir_name: str) -> list[str]: 
    """
    Takes a directory and returns a list with all the names in it in sorted order
    """
    file_names = sorted([os.path.splitext(name)[0] for name in os.listdir(dir_name)], key=lambda x: (int(x.split('_')[0]), x.split('_')[1], x.split('_')[2]))
    return file_names

#####################

def check_file_lists(file_list1: list[str], file_list2: list[str]) -> None: 
    """
    Takes to list of files and checks that the lists contains the same elements in the same order
    """
    assert len(file_list1) == len(file_list2)

    for n in range(len(file_list2)): 
        try: 
            assert file_list1[n] == file_list2[n]
        except: 
            raise ValueError(f"ERROR: {file_list1[n]} is not equal to {file_list2[n]}")

def check_content(file_names: list[str], sequences_path: str, structures_path: str) -> None: 
    """
    Files in file_names should be named as length_type_rest
    For every file in file_names in checks that the length of the sequence and structure correspons to the name
    Assumes that the paths contains the same files and the files in file_names
      
    Args: 
        - files_names: a list of the files that should be in the directories
        - sequence_path: name of path that contains the files with sequences
        - structures_path: name of path that contains the files with structures
    """
    for n in range(len(file_names)): 
        file_name = file_names[n]
        length = int(file_name.split("_")[0])
        sequence, _ = read_fasta(os.path.join(sequences_path, file_name + ".fasta"))
        structure, _ = read_dbn_file(os.path.join(structures_path, file_name + ".dbn"))

        if (length != len(sequence)) & (length != len(structure)):
            raise ValueError(f"ERROR: {file_name} named wrong!")
        elif len(sequence) != len(structure): 
            raise ValueError(f"ERROR: Structure and sequence for {file_name} have different lengths")

def check_files(sequence_path: str, structures_path: str) -> None:
    """
    Collects the checks of files 
    Takes a path containing files with sequences and one containing files with structures. 
    Checks that the direcotries contains the same files and the content of the files are constistent with the name in regards to length
    """    
    filelist1 = get_path_list(sequence_path)
    filelist2 = get_path_list(structures_path)

    print(f"Cheking {sequence_path} and {structures_path}", file=sys.stderr)

    check_file_lists(filelist1, filelist2)

    check_content(filelist1, sequence_path, structures_path)



def main() -> None:
    """
    Is called as check_files *structure_paths
    It checks the consistency between the folder sequences and all provided structure paths in the folder structure
    """ 
    arguments = sys.argv[1:]
    sequence_path = "../sequences"
    for arg in arguments:
        structures_path = "../structures/" + str(arg)

        check_files(sequence_path, structures_path)


if __name__ == "__main__": 
    main()