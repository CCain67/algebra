import pytest

import galois

from group_theory.common_groups.residue_class_groups import *
from group_theory.common_groups.permutation_groups import *
from group_theory.common_groups.matrix_groups import *
from group_theory.common_groups.misc_groups import *

def test_cyclic_group_order():
    assert CyclicGroup.as_residue_class_group(31).order == 31

def test_multiplicative_residue_class_group_order():
    assert GroupOfUnits.as_residue_class_group(16).order == galois.euler_phi(16)

def test_symmetric_group_order():
    assert SymmetricGroup.as_permutation_group(4).order == 24
    assert SymmetricGroup.as_matrix_group(4).order == 24

def test_dihedral_group():
    assert DihedralGroup.as_permutation_group(7).order == 14

def get_gl_n_p_order(n,q):
    k = 0
    prod = 1 
    while k<=n-1:
        prod *= (q**n-q**k)
        k += 1
    return prod

def test_gl_n_p_order():
    assert GeneralLinear.over_finite_field(2,3,1).order == get_gl_n_p_order(2,3)

def test_sl_n_p_order():
    assert SpecialLinear.over_finite_field(2,3,1).order == get_gl_n_p_order(2,3)/2

def test_heisenberg_order():
    assert HeisenbergGroup.over_finite_field(3,5,1).order == 5**3