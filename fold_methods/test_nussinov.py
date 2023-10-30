import pytest
import numpy as np
from io import StringIO

from nussinov_expanded import read_fasta, prepare_input, pairing_score, bifurcating_score, fill_S, find_optimal, backtrack, fold_RNA, db_to_file, read_parameter

def test_read_fasta(): 
    """
    Test that read_fasta return the sequence as a string. 
    If the fasta file contains more than one sequence a ValueError is raised.
    """
    assert read_fasta(StringIO('>seq1\naaaccctcg')) == ('aaaccctcg', 'seq1')

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

def test_pairing_score(): 
    """
    """
    parameters = read_parameter("pairing_parameters.csv")
    S = np.zeros((20, 20))
    
    sequence1 = "CCCCCCCCCCGGGGGGGGGG"
    sequence2 = "CCCCCCCCCCCGCGCGCGCG"

    for i in range(9, -1, -1):  
        j = 19-i
        if i % 2 == 0:
            assert pairing_score(i, j, S, sequence2, parameters) == 0
        else: 
            assert pairing_score(i, j, S, sequence2, parameters) == float('inf')

        assert pairing_score(i, j, S, sequence1, parameters) == parameters.at[sequence1[i]+sequence1[j], sequence1[i+1]+sequence1[j-1]]

def test_bifurcating_score(): 
    """
    Tests that the correct value is found for the bifurcating loops
    """
    S = np.array([[0]*10, [1]*10, [2]*10, [3]*10, [4]*10, [5]*10, [6]*10, [7]*10, [8]*10, [9]*10])
    
    i_list = [i for i in range(0, 4)]
    j_list = [j for j in range(9, 6, -1)]

    energies = [2, 4, 6]
    ij = [1, 2, 3]
    
    for n in range(3): 
        i, j = i_list[n], j_list[n]
        assert bifurcating_score(i, j, S) == (energies[n], ij[n])

def test_fill_S(): 
    parameters = read_parameter("pairing_parameters.csv")
    S_true = np.array([[0, 0, 0, 0, 0, 0, 0, -1.3, -4.2],
                       [0, 0, 0, 0, 0, 0, 0, -1.3, -2.9],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                       [0, 0, 0, 0, 0, 0, 0, 0, 0]])
    sequence = "GGGAAAUCC"

    assert fill_S(sequence, parameters).all() == S_true.all()

def test_find_optimal(): 
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

def test_fold_RNA(): 
    parameters = read_parameter("pairing_parameters.csv")
    S = np.array([[0, 0, 0, 0, 0, 0, 0, -1.3, -4.2],
                       [0, 0, 0, 0, 0, 0, 0, -1.3, -2.9],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0, 0, 0, 0, 0], 
                       [0, 0, 0, 0, 0, 0, 0, 0, 0]])
    sequence = "GGGAAAUCC"
    
    assert fold_RNA(S, sequence, parameters) == '(((...)))' 