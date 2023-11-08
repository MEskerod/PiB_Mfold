from gwf import Workflow
from gwf import AnonymousTarget

import os

gwf = Workflow()

def analysis(): 
    inputs = []
    outputs = [os.path.join("results", "time_table.csv"),
               os.path.join("results", "synthethic_times.csv"),
               os.path.join("results", "distance_table.csv"),
               os.path.join("results", "distances_newest.jpeg"),
               os.path.join("results", "distances_true.jpeg"),
               os.path.join("results", "Nussinov_times.jpeg"),
               os.path.join("results", "original_vs_new_times.jpeg"),
               os.path.join("results", "synthetic_vs_real_times.jpeg")]
    options = {"memory":"8gb", "walltime":"32:00:00"}
    spec = """cd analysis
    
    python3 analysis.py""".format()
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

gwf.target_from_template(name = "Mfold_analysis", template = analysis())