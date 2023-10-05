import subprocess

### FUNCTIONS FOR DIFFERENT PROGRAMS ###
"""
The functions have to take file and output directory as parameters. 
Each function runs the functionality of the program on 'file' and outputs the resulting dot bracket structure to a .txt file
"""

def run_Mfold_web(file, outdir):
    command = f"python3 cmd_Mfold_db.py {file} {outdir}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()

def run_Mfold_orginal(file, outdir):
    outputfile = os.path.join(outdir, file)
    command = f"python3 Mfold/main.py -f {file} -lp 1988 -o {outputfile}"
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()