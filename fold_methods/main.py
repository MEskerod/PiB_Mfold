### FILES USED BY MULTIPLE FILES IN FOLDER ###  

def db_to_file(sequence, db, filename, name): 
    with open(filename, 'w') as f: 
        f.write(f"#Name: {name}\n")
        f.write(f"#Length: {len(sequence)}\n")
        f.write(sequence + "\n")
        f.write(db + "\n")