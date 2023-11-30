import subprocess, os, sys

### FUNCTIONS FOR DIFFERENT PROGRAMS ###
"""
The functions have to take file and output directory as parameters. 
Each function runs the functionality of the program on 'file' and outputs the resulting dot bracket structure to a .txt file
"""

def run_Mfold_web(file: str, outdir: str) -> None:
    """
    Runs the version of Mfold downloaded from unafold website

    Args: 
        - file: fasta file containing sequence that should be folded
        - outdir: directory, where structure should be outputted to
    """
    command = f"python3 ../fold_methods/cmd_Mfold_db.py {file} {outdir}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    if p.returncode != 0:
        print(f"Command failed for Mfold Web for {file} with return code {p.returncode}", file = sys.stderr)
        print("Error output:", file = sys.stderr)
        print(p.stderr.read().decode('utf-8'), file = sys.stderr)


def run_Mfold_orginal(file: str, outdir: str) -> None:
    """
    Runs the orginal version of Mfold that was implemeted

    Args: 
        - file: fasta file containing sequence that should be folded
        - outdir: directory, where structure should be outputted to
    """
    outputfile = os.path.join(outdir, os.path.splitext(os.path.basename(file))[0])
    command = f"python3 ../Mfold/main.py -f {file} -lp 1988 -o {outputfile}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    if p.returncode != 0:
        print(f"Command failed for Mfold original for {file} with return code {p.returncode}", file = sys.stderr)
        print("Error output:", file = sys.stderr)
        print(p.stderr.read().decode('utf-8'), file = sys.stderr)

def run_Mfold_newest(file: str, outdir: str) -> None: 
    """
    Runs the newest version of Mfold that was implemented

    Args: 
        - file: fasta file containing sequence that should be folded
        - outdir: directory, where structure should be outputted to
    """
    outputfile = os.path.join(outdir, os.path.splitext(os.path.basename(file))[0])
    command = f"python3 ../Mfold/main.py -f {file} -b -a -c -o {outputfile}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    if p.returncode != 0:
        print(f"Command failed for newest version of Mfold for {file} with return code {p.returncode}", file = sys.stderr)
        print("Error output:", file = sys.stderr)
        print(p.stderr.read().decode('utf-8'), file = sys.stderr)

def run_Nussinov(file: str, outdir: str) -> None: 
    """
    Runs the implemented version of Nussinov, where energy parameters for stacking is used

    Args: 
        - file: fasta file containing sequence that should be folded
        - outdir: directory, where structure should be outputted to
    """
    outputfile = os.path.join(outdir, os.path.splitext(os.path.basename(file))[0])
    command = f"python3 ../fold_methods/nussinov_expanded.py -f {file} -o {outputfile}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

    if p.returncode != 0:
        print(f"Command failed for Nussinov for {file} with return code {p.returncode}", file = sys.stderr)
        print("Error output:", file = sys.stderr)
        print(p.stderr.read().decode('utf-8'), file = sys.stderr)