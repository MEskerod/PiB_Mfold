import pytest
import numpy as np
from io import StringIO

from Mfold import(
    read_fasta,
    prepare_input,
    read_parameters,
    loop_greater_10,
    find_E1,
    stacking,
    bulge_loop,
    interior_loop,
    find_E2,
    find_E3,
    penta_nucleotides,
    compute_V,
    find_E4,
    compute_W,
    fold_rna,
    find_optimal,
    backtrack
)

def test_read_fasta(): 
    """
    Test that read_fasta return the sequence as a string. 
    If the fasta file contains more than one sequence a ValueError is raised.
    """
    assert read_fasta(StringIO('>seq1\naaaccctcg')) == 'aaaccctcg'

    with pytest.raises(ValueError): 
        read_fasta(StringIO('>seq1\nAATTAAT\n>seq2\naattaat\n>seq3\nAaTtAaT'))

def test_input(): 
    """
    Test that prepare input can: 
    - handle upper and lower case
    - change T to U 
    - raise an error if non-valid letters are present in the sequence
    """
    assert prepare_input('AACGC') == 'AACGC'
    assert prepare_input('aacgc') == 'AACGC'
    assert prepare_input('aacgc\n') == 'AACGC'
    assert prepare_input('ttt') == 'UUU'

    with pytest.raises(ValueError): 
        prepare_input('aacgy')

def test_parameters(): 
    """
    """
    assert 2==2

def test_loop_greater_10(): 
    """
    Tests that the correct value is calculated for loops of sizes greater than 10
    """
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    
    IL = [6.35, 6.76, 7.05, 7.28, 7.47]
    BL =[5.65, 6.05, 6.35, 6.58, 6.77]
    HL = [6.02, 6.42, 6.71, 6.94, 7.13]

    sizes = [11, 16, 21, 26, 31]

    for n in range(5): 
        assert round(loop_greater_10("IL", sizes[n], loop),2) == IL[n]
        assert round(loop_greater_10("BL", sizes[n], loop),2) == BL[n]
        assert round(loop_greater_10("HL", sizes[n], loop),2) == HL[n]

def test_E1(): 
    """
    The energy of a hairpin loop should corespond to the size of the loop 
    If size >10 the energy should be caculated
    """
    parameters = [4.5,5.5,4.9,5.1,5.2,5.5,5.8,5.9]
    
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    i = 0 
    for n in range(8): 
        j = n+4
        assert find_E1(i, j, loop) == parameters[n] 

    for j in range(12,32,5): 
        assert find_E1(i, j, loop) == loop_greater_10("HL", j-i-1, loop)

def test_stacking(): 
    """
    Tests that looking up base pairs in the table returns the correct value
    Furthermore tests that the value added is V[i+1, j-1]
    """
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    
    V = np.array([[1,1,1,1],[1,1,0,1], [1,1,1,1], [1,1,1,1]])


    sequences05 = ['AGUU', 'UUGA', 'GAUU', 'GGUU', 'GUGU', 'UUAG', 'UGUG', 'UUGG']
    sequences07 = ['AUGU', 'UGUA', 'GUAU', 'UAUG']
    sequences09 = ['AAUU', 'UUAA', 'AUAU']
    sequences11 = ['UAUA']
    sequences13 = ['GGUC', 'UCGG']
    sequences15 = ['CGUG', 'CUGG', 'GGCU', 'UGCG']
    sequences17 = ['AGCU', 'CUAG']
    sequences18 = ['UGCA', 'CAUG']
    sequences19 = ['GUGC', 'GCGU']
    sequences20 = ['CGCG']
    sequences21 = ['ACGU', 'GUAC']
    sequences23 = ['UCGA', 'GAUC']
    sequences29 = ['CCGG', 'GGCC']
    sequences34 = ['GCGC']

    scores = [-0.5, -0.7, -0.9, -1.1, -1.3, -1.5, -1.7, -1.8, -1.9, -2.0, -2.1, -2.3, -2.9, -3.4]
    pairs = [sequences05, sequences07, sequences09, sequences11, sequences13, sequences15, sequences17, 
             sequences18, sequences19, sequences20, sequences21, sequences23, sequences29, sequences34]

    i, j = 0, 3

    for n in range(14): 
        for sequence in pairs[n]: 
           assert round(stacking(i, j, V, stack, sequence),2) == scores[n] 

def test_bulgeloop(): 
    assert 2==2

def test_interiorloop(): 
    assert 2==2

def test_bifurcation(): 
    assert 2==2

def test_computeV(): 
    assert 2==2

def test_computeW(): 
    assert 2==2

def test_foldRNA(): 
    assert 2==2

def test_optimal(): 
    assert 2==2

def test_backtrack(): 
    assert 2==2

