import pytest, sys, subprocess
import numpy as np
from io import StringIO

from fold_functions import (make_asymmetric_penalty, loop_greater_10, stacking, bulge_loop_3end, bulge_loop_5end, interior_loop, 
                            find_E1, find_E3, find_E4, compute_V, fold_rna, find_optimal, backtrack)

from help_functions import (read_fasta, prepare_input, read_general_parameters, parse_asymmetry_parameters)

### TEST HELPER FUNCTIONS ###
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

def test_parse_asymmetry_parameters(): 
    assert parse_asymmetry_parameters("[0.4, 0.3, 0.2, 0.1]; 3") == ([0.4, 0.3, 0.2, 0.1], 3)
    assert parse_asymmetry_parameters("[0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3]; 123") == ([0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3], 123) 

    with pytest.raises(Exception): 
        parse_asymmetry_parameters("[0.4, 0.3, 0.2, 0.1]")
    
    with pytest.raises(Exception): 
        parse_asymmetry_parameters("0.4, 0.3, 0.2, 0.1; 3")
    
    with pytest.raises(Exception): 
        parse_asymmetry_parameters("[0.4, 0.3, 0.2, 0.1], 3")
    
    with pytest.raises(Exception): 
        parse_asymmetry_parameters("[0.4, 0.3, 0.2, 0.1] 3")
    
    with pytest.raises(Exception): 
        parse_asymmetry_parameters("[0.4, 0.3, 0.2, 0.1]; 3 4")


### TEST FOLDING FUNCTIONS ###
def test_loop_greater_10(): 
    """
    Tests that the correct value is calculated for loops of sizes greater than 10
    """
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    
    IL = [6.4, 6.81, 7.1, 7.33, 7.52]
    BL =[5.6, 6.01, 6.3 , 6.53, 6.72]
    HL = [6., 6.41, 6.7, 6.93, 7.12]

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
    
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    i = 0 
    for n in range(8): 
        j = n+4
        assert find_E1(i, j, loop) == parameters[n] 

    for j in range(12,32,5): 
        assert find_E1(i, j, loop) == round(loop_greater_10("HL", j-i-1, loop), 5)

def test_stacking(): 
    """
    Tests that looking up base pairs in the table returns the correct value
    Furthermore tests that the value added is V[i+1, j-1]
    """
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    
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
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    
    V = np.full((34,34), float('inf'))
    V[11, :] = 0

    sequence = 'CCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGG'

    i_list, j_list = [x for x in range(8)], [x for x in range(33, 25, -1)]
    
    for n in range(8): 
        i, j = i_list[n], j_list[n]
        assert bulge_loop_5end(i, j, V, loop, stack, sequence, stacking=True) == (loop.at[10-i, "BL"], 11)

def test_bulgelopp_j():
    """
    Tests that bulge loops on the 3' end of size 2-10 returns the correct value
    """
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    
    V = np.full((34,34), float('inf'))
    V[:, 22] = 0

    sequence = 'AAAAAAAAAAAAAAAAAUUUUUUUUUUUUUUUUU'

    i_list, j_list = [x for x in range(8)], [x for x in range(33, 25, -1)]
    
    for n in range(8): 
        i, j = i_list[n], j_list[n]
        assert bulge_loop_3end(i, j, V, loop, stack, sequence, stacking=True) == (loop.at[10-i, "BL"], 22)

def test_bulgeloop_size1(): 
    """
    Tests that bulge loops of size 1 returns the loop parameter, stacking parameter and V[i-2, j-1]/V[i-1, j-2]
    """
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")

    V = np.full((10,10), float('inf'))
    V[2, 8], V[2, 6] = 0, 0

    sequences = ['AAAAAUUUUU', 'CCCCCGGGGG', 'UUUUUGGGGG']

    #Bulge on i
    i, j = 0, 9

    for sequence in sequences: 
        bp = sequence[i] + sequence[j]
        assert bulge_loop_5end(i, j, V, loop, stack, sequence, stacking=True) == (loop.at[1, "BL"] + stack.at[bp, bp], i+2)
    
    #Bulge on j
    i, j = 1, 8

    for sequence in sequences: 
        bp = sequence[i] + sequence[j]
        assert bulge_loop_3end(i, j, V, loop, stack, sequence, stacking=True) == (loop.at[1, "BL"] + stack.at[bp, bp], j-2)

def test_interiorloop_symmetric(): 
    """
    Test that the energy of symmetric interior loops is the same as given in the table
    """
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    asymetric_penalty_function = make_asymmetric_penalty([0.4, 0.3, 0.2, 0.1], 3)

    V = np.full((34,34), float('inf'))
    V[6,27] = 0

    sequence = 'CCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGG'

    i_list = [x for x in range(5)]
    j_ist = [x for x in range(33, 28, -1)]

    for n in range(5): 
        assert interior_loop(i_list[n], j_ist[n], V, loop, sequence, asymetric_penalty_function, False, True) == (loop.at[10-(n*2), "IL"], (6, 27))

def test_interiorloop_asymmetric(): 
    """
    Test that the energy of asymmetric interior loops is the same as given in the table + the penalty
    """
    loop, stack = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    asymetric_penalty_function = make_asymmetric_penalty([0.4, 0.3, 0.2, 0.1], 3)

    sequence = 'CCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGG'

    penalties = [0.1*1, 0.2*2, 0.3*3, 0.4*4]

    #5' loop > 3' end loop
    i = 0
    j_ist = [x for x in range(33, 29, -1)]
    V = np.full((34,34), float('inf'))
    V[6,28] = 0
    for n in range(4): 
        assert interior_loop(i, j_ist[n], V, loop, sequence, asymetric_penalty_function, False, True) == (round(loop.at[10-n-1, "IL"] + penalties[n], 5), (6, 28))
    
    #3' loop > 5' end loop
    i_list = [x for x in range(4)]
    j = 33
    V = np.full((34,34), float('inf'))
    V[5,27] = 0
    for n in range(4): 
        assert interior_loop(i_list[n], j, V, loop, sequence, asymetric_penalty_function, False, True) == (round(loop.at[10-n-1, "IL"] + penalties[n], 5), (5, 27))    

def test_bifurcation(): 
    """
    Tests that the correct value is found for the bifurcating loops
    """
    W = np.array([[0]*10, [1]*10, [2]*10, [3]*10, [4]*10, [5]*10, [6]*10, [7]*10, [8]*10, [9]*10])
    
    i_list = [i for i in range(0, 4)]
    j_list = [j for j in range(9, 6, -1)]

    energies = [4, 6, 8]
    ij = [(2,3), (3,4), (4,5)]
    
    for n in range(3): 
        i, j = i_list[n], j_list[n]
        assert find_E3(i, j, W) == (energies[n], ij[n])

def test_E4():
    """
    Tests that the corret value is found for E4, which is i and j both base pair, but with different bases
    """
    
    W = np.array([[0]*10, [1]*10, [2]*10, [3]*10, [4]*10, [5]*10, [6]*10, [7]*10, [8]*10, [9]*10])
    
    i_list = [i for i in range(0, 4)]
    j_list = [j for j in range(9, 5, -1)]

    energies = [2, 4, 6, 8]
    ij = [(1,2), (2,3), (3,4), (4,5)]
    
    for n in range(4): 
        i, j = i_list[n], j_list[n]
        assert find_E4(i, j, W) == (energies[n], ij[n])

def test_computeV(): 
    """
    Tests that if i and j base pair a value is calculated for V[i,j] and otherwise the value should be inf
    """
    basepairs = ['AU', 'UA', 'CG', 'GC', 'GU', 'UG']
    parameters = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    asymetric_penalty_function = make_asymmetric_penalty([0.4, 0.3, 0.2, 0.1], 3)
    sequence = "ACGUACGUACGUACGUACGUACGU"
    N = len(sequence)
    V, W = np.zeros((N, N)), np.zeros((N, N))

    for l in range(5, N): #Computes the best score for all subsequences that are longer than 5 nucleotides
        for i in range(0, N-5): 
            j = i+l
            if j < N: 
                compute_V(i, j, W, V, sequence, parameters, asymetric_penalty_function, True, False, True)
                if sequence[i]+sequence[j] in basepairs:
                    assert V[i,j] != float('inf')
                else: 
                    assert V[i,j] == float('inf')

def test_foldRNA(): 
    """
    Tests that the correct matrices are outputted for one sequence
    """
    parameters = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    asymetric_penalty_function = make_asymmetric_penalty([0.4, 0.3, 0.2, 0.1], 3)

    W, V = np.full((12, 12), float('inf')), np.full((12, 12), float('inf'))
    W[0, 5:] = [4.5, 3.6, 3.6, 3.6, 2.5, 1.2, 0.5]
    W[1, 5:] = [4.5, 4.5, 4.5, 4.5, 2.5, 1.2, 0.5]
    W[2, 8:] = [4.5, 2.5, 1.2, 1.2]
    W[3, 8:] = [4.5, 2.5, 2.5, 2.5]
    W[4, 8:] = [4.5, 4.5, 4.5, 4.5]
    W[5, 9:] = [4.5, 4.5, 4.5]
    W[6, 10:] = [4.5, 4.5]

    V[0, 5:7] = [5.5, 3.6]
    V[0, 11] = 4.4
    V[1, 5:7] = [4.5, 5.5]
    V[1, 11] = 0.5
    V[2, 9:11] = [5.1, 1.2]
    V[3, 9:11] = [2.5, 5.1]
    V[4, 8] = 4.5 
    V[4, 11] = 5.0
    V[5, 9:11] = [4.5, 5.5]
    V[6, 10] = 4.5
    
    true_W, true_V = fold_rna("AAUCGUUCCGGU", parameters, asymetric_penalty_function, True, False, True)
    assert true_W.all() == W.all() 
    assert true_V.all() == V.all()

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
    """
    Tests that backtrack is done correctly trough for one sequence
    """
    parameters = read_general_parameters("parameters/loop_1989.csv", "parameters/stacking_1988.csv")
    assymetric_penalty_function = make_asymmetric_penalty([0.4, 0.3, 0.2, 0.1], 3)
    W, V = fold_rna("AAUCGUUCCGGU", parameters, assymetric_penalty_function, True, False, True)

    assert backtrack(W, V, parameters, "AAUCGUUCCGGU", assymetric_penalty_function, True, False, True) == '.((((...))))'

def test_argparser():
    """
    Tests that the arguments parsed are handled appropiatly
    """ 

    #Other parameters
    result = subprocess.Popen([sys.executable, "main.py", "-i", "AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU", "-lp", "1988"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = result.communicate()
    assert "Energy of optimal fold is -56.3 kcal/mol" in stdout.decode('utf-8')
    assert "..(((((((((....((((((((..(((..((((..((.....))..))))..)))...))))))))..((((((...((.((((((((.....))))))))..)).))))))))))))))).." in stdout.decode('utf-8')
    assert '' == stderr.decode('utf-8')

    #With bulge stacking
    result = subprocess.Popen([sys.executable, "main.py", "-i", "AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU", "-b"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = result.communicate()
    assert "Energy of optimal fold is -50.0 kcal/mol" in stdout.decode('utf-8')
    assert "..(((((((((....((((((((........((((.....))))(((......)))...))))))))..(((((((..((.((((((((.....))))))))..)))))))))))))))))).." in stdout.decode('utf-8') 
    assert '' == stderr.decode('utf-8')
    
    #Asymmetric
    result = subprocess.Popen([sys.executable, "main.py", "-i", "AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU", "-a", "-b"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = result.communicate()
    assert "Energy of optimal fold is -49.6 kcal/mol" in stdout.decode('utf-8')
    assert "..(((((((((....((((((((........((((.....))))(((......)))...))))))))..(((((((..((.((((((((.....))))))))..)))))))))))))))))).." in stdout.decode('utf-8')
    assert '' == stderr.decode('utf-8')

    #Other assymetry parameters
    result = subprocess.Popen([sys.executable, "main.py", "-i", "AUGCCGACGGUCAUAGGACGGGGGAAACACCCGGACUCAUUCCGAACCCGGAAGUUAAGCCCCGUUCCGUCCCGCACAGUACUGUGUUCCGAGAGGGCACGGGAACUGCGGGAACCGUCGGCUU", "-a", "-A", "[0.7, 0.6, 0.4, 0.2, 0.1]; 6", "-b"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = result.communicate()
    assert '' == stderr.decode('utf-8')
    assert "Energy of optimal fold is -49.5 kcal/mol" in stdout.decode('utf-8')
    assert "..(((((((((....((((((((........((((.....))))(((......)))...))))))))..(((((((......(((((((.....)))))))(...))))))))))))))))).." in stdout.decode('utf-8') 