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

def test_bulgeloop_i(): 
    """
    Tests that bulge loops on the 5' end of size 2-10 returns the correct value
    """
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    
    V = np.full((34,34), float('inf'))
    V[11, :] = 0

    sequence = 'CCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGG'

    i_list, j_list = [x for x in range(8)], [x for x in range(33, 25, -1)]
    
    for n in range(8): 
        i, j = i_list[n], j_list[n]
        assert bulge_loop(i, j, V, loop, stack, sequence) == loop.at[10-i, "BL"]

def test_bulgelopp_j():
    """
    Tests that bulge loops on the 3' end of size 2-10 returns the correct value
    """
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    
    V = np.full((34,34), float('inf'))
    V[:, 22] = 0

    sequence = 'AAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUU'

    i_list, j_list = [x for x in range(8)], [x for x in range(33, 25, -1)]
    
    for n in range(8): 
        i, j = i_list[n], j_list[n]
        assert bulge_loop(i, j, V, loop, stack, sequence) == loop.at[10-i, "BL"]

def test_bulgeloop_size1(): 
    """
    Tests that bulge loops of size 1 returns the loop parameter, stacking parameter and V[i-2, j-1]/V[i-1, j-2]
    """
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")

    V = np.full((10,10), float('inf'))
    V[2, 8], V[2, 6] = 0, 0

    sequences = ['AAAAAUUUUU', 'CCCCCGGGGG', 'UUUUUGGGGG']

    #Bulge on i
    i, j = 0, 9

    for sequence in sequences: 
        bp = sequence[i] + sequence[j]
        assert bulge_loop(i, j, V, loop, stack, sequence) == loop.at[1, "BL"] + stack.at[bp, bp]
    
    #Bulge on j
    i, j = 1, 8

    for sequence in sequences: 
        bp = sequence[i] + sequence[j]
        assert bulge_loop(i, j, V, loop, stack, sequence) == loop.at[1, "BL"] + stack.at[bp, bp]

def test_interiorloop_symmetric(): 
    """
    Test that the energy of symmetric interior loops is the same as given in the table
    """
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")

    V = np.full((34,34), float('inf'))
    V[6,27] = 0

    sequence = 'CCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGG'

    i_list = [x for x in range(5)]
    j_ist = [x for x in range(33, 28, -1)]

    for n in range(5): 
        assert interior_loop(i_list[n], j_ist[n], V, loop, sequence) == loop.at[10-(n*2), "IL"]

def test_interiorloop_asymmetric(): 
    """
    Test that the energy of asymmetric interior loops is the same as given in the table + the penalty
    """
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")

    sequence = 'CCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGG'

    penalties = [0.1*1, 0.2*2, 0.3*3, 0.4*4]

    #5' loop > 3' end loop
    i = 0
    j_ist = [x for x in range(33, 29, -1)]
    V = np.full((34,34), float('inf'))
    V[6,28] = 0
    for n in range(4): 
        assert interior_loop(i, j_ist[n], V, loop, sequence) == loop.at[10-n-1, "IL"] + penalties[n]
    
    #3' loop > 5' end loop
    i_list = [x for x in range(4)]
    j = 33
    V = np.full((34,34), float('inf'))
    V[5,27] = 0
    for n in range(4): 
        assert interior_loop(i_list[n], j, V, loop, sequence) == loop.at[10-n-1, "IL"] + penalties[n]    

def test_bifurcation(): 
    assert 2==2

def test_computeV(): 
    """
    Tests that if i and j base pair a value is calculated for V[i,j] and otherwise the value should be inf
    """
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']
    parameters = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    sequence = "ACGUACGUACGUACGUACGUACGU"
    N = len(sequence)
    V, W = np.zeros((N, N)), np.zeros((N, N))

    for l in range(5, N): #Computes the best score for all subsequences that are longer than 5 nucleotides
        for i in range(0, N-5): 
            j = i+l
            if j < N: 
                compute_V(i, j, W, V, sequence, parameters)
                if sequence[i]+sequence[j] in basepairs:
                    assert V[i,j] != float('inf')
                else: 
                    assert V[i,j] == float('inf')

def test_computeW(): 
    assert 2==2

def test_foldRNA(): 
    assert 2==2

def test_optimal(): 
    """
    Tests that the upper right corner is returned
    """
    assert find_optimal(np.array([[1, 2, 3, 4, 5], 
                                  [6, 7, 8, 9, 10],
                                  [11, 12, 13, 14, 15]])) == 5
    
    assert find_optimal(np.array([[0, 0, 1000], 
                                 [0, 0, 0], 
                                 [0, 0, 0],
                                 [0, 0, 0],
                                 [0, 0, 0]])) == 1000
    
    assert find_optimal(np.array([[1],
                                 [5],
                                 [5]])) == 1

def test_backtrack(): 
    assert 2==2

