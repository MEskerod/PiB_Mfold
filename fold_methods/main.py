import pandas as pd

### FILES USED BY MULTIPLE FILES IN FOLDER ###  

def db_to_file(sequence, db, filename, name): 
    with open(filename, 'w') as f: 
        f.write(f"#Name: {name}\n")
        f.write(f"#Length: {len(sequence)}\n")
        f.write(sequence + "\n")
        f.write(db + "\n")


def read_parameter(file_stacking): 
    """

    """
    try:
        stacking = pd.read_csv(file_stacking, index_col=0)
    except FileNotFoundError:
        raise FileNotFoundError("Parameter file not found")
    except pd.errors.EmptyDataError:
        raise ValueError("Parameter file is empty or in an unexpected file format")
    return stacking