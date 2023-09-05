import pytest

from Mfold import(
    read_fasta,
    prepare_input,
    read_parameters,
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


def test_E1(): 
    """
    The energy of a hairpin loop should corespond to the size of the loop 
    If size >? do ?
    """
    parameters = [4.5,5.5,4.9,5.1,5.2,5.5,5.8,5.9]
    
    loop, stack = read_parameters("loop_improved.csv", "pairing_parameters.csv")
    i = 0 
    for n in range(8): 
        j = n+4
        assert find_E1(i, j, loop) == parameters[n] 

def test_stacking(): 
    assert 2==2 

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

