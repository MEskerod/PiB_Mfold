# Mfold implemented in Python 
Mfold algorithm for RNA folding implemented as described in:   
** Add articles ** 

## Running from command line

```
python3 main (???)   -f 'filename.fasta'/-i 'sequence' -o 'outputfile'
```

**Input sequence**   
The input sequence can be either a file or inputted in the command.    
If input is in a file it should be a fasta file with a single sequence and supplied as -f 'filename.fasta'.    
Else the sequence should be inputted as -i 'sequence' and the fold will be named 'user inputtet sequence' in the output. 

**Output**   
Default is to output the dot bracket structure to the terminal.    
If the structure should be outputted to a file this can be specified using -o followed by the output file name without extension. 

The fold is outputted in the format of a .dbn file.    
```
#Name: 'sequence name'
#Length: 'sequence length'
SEQUENCE 
DOT BRACKET STRUCTURE
```



### Optional arguments
Default parameters will give the newest of the implementedd versions. 

**PARAMETERS**
```
-lp --loop_parameters 'year' (default = 1989)
```
Makes it possible to chose from different predefined sets of loops parameters.    
The choices are: 
- 1988 (described in __)
- 1989 (described in __)

**PENALTY FOR ASYMMETRIC INTERIOR LOOPS**   
Asymmetric interior loops are not as stable as symmetric and i some versions of Mfold this is penalized according to the size of the loop. 

```
-a --asymmetric
```
If -a/--asymmetric is given a penalty for asymmetric loops is added. A default parameters for the penalty function is used. It can be specified using the flag -A/--asymmetry_parameters (see blow)

**CONTINUED BASE PAIR STACKING FOR BULGE LOOPS OF SIZE 1**   
In the Turner article from 1989 continued stacking for bulge loops of size 1 was described. When the bulge is not larger than 1 the favorable energy of stacking is retained and the parameter for base pairing is added to the one for the loop. 
```
-b --bulge_stacking 
```
If -b/--bulge_stacking is given stacking is retained. 


**PENALTY FOR GU AND AU BASE PAIRS AS CLOSING BASE PAIRS FOR INTERIOR LOOPS**   
Interior loops closed by one or two AU or GU base pairs is not as stable as loops enclosed by GC base pairs. Interior loops closed by one or two AU/GU base pairs are penaized by additional 0.9 kcal/mol and 1.8 kcal/mol respectively 
```
-c --closing_penalty
```
If -c/--closing_penalty is given penalty is added for closing GU/AU base pairs.

**PARAMETERS FOR PENALTY FUNCTION FOR ASYMMETRIC INTERIOR LOOPS**   
The penalty that is added for asymmetric loops is calculated as the minimum of a penalty maximum and N*f(M), where 
- N is |N1-N2| and N1 and N2 is the size of the loop on each strand
- f is a function of the possible M. Here it is supplied as a list
- M is the minimum of the length of f, N1 and N2   

The parameters for the asymmetric penalty function is supplied as: 
```
-A --asymmetry_parameters "f; penalty_max" (default = "[0.4, 0.3, 0.2, 0.1]; 3")
```
**NOTE:** The asymmetry parameters must be supplied enclosed in " " and with *f* (list) seperated from *penalty_max* (int) by a semicolon

**SUPPLY DIFFERENT STACKING OR LOOP PARAMETERS**
```
-lf --loop_file 'file_name.csv'

-sf --stacking_file 'file_name.csv'
```

If different parameter files are needed they can be supplied using -lf or -sf.    
The parameter files should be .csv files. 

The file with loop parameters has to have a header of *IL,BL,HL* refering to interior loop, bulge loop and hairpin loop. It should contain 11 rows with the energy for sizes from 0-10. If the loop size is not allowed (ex. HL, size 3) it is left empty.   

The file containing stacking parameters

## Arguments for versions described in articles 
**Original version with parameters from Turner 1988:**   

```
python3 main -f/-i 'input' -lp 1988 (-o 'output_file')
```

**Version described by Turner in 1988:**   
```
python3 main -f/-i 'input' -lp 1988 -b -a -c -A "[0.7, 0.6, 0.4, 0.2, 0.1]; 6" (-o 'output_file')
```

**Version described by Jaeger and Zuker in 1989 without dangling ends and pentanucleotides:**  
```
python3 main -f/-i 'input' -b -a -c (-o 'output_file')
``` 
