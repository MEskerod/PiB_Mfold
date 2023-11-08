import os, sys

from general import make_dir, get_path_list
from running_fold import run_Mfold_web

def main(): 

    print(f"\nFOLDING WITH MFOLD WEB", file=sys.stderr)
    
    output_dir = make_dir(os.path.join("..", "structures", "MfoldWebversion" + "_structures")) 

    file_list = get_path_list("../sequences")

    for index, file in enumerate(file_list): 
        print(f"Fold {file} - {index + 1}/{len(file_list)}", file=sys.stderr)
        run_Mfold_web(file, output_dir)

if __name__ == '__main__': 
    main()