import pytest 

from group_theory.base.groups import Subgroup

from group_theory.common_groups.permutation_groups import SymmetricGroup
from group_theory.common_groups.misc_groups import DihedralGroup

def test_center_dihedral():
    assert DihedralGroup.as_permutation_group(7).center().order == 1
    assert DihedralGroup.as_permutation_group(8).center().order == 2

def test_normalizer():
    G = SymmetricGroup.as_permutation_group(5)
    S = G.subgroup_generated_by([G.canonical_generators[0]])
    N = G.normalizer(S)
    C = G.centralizer(S)
    S = Subgroup(S.elements,N)
    assert S.is_normal
    assert C<N


def test_derived_series():
    S_4 = SymmetricGroup.as_permutation_group(4)
    S_5 = SymmetricGroup.as_permutation_group(5)
    assert S_4.is_solvable()
    assert not S_5.is_solvable()

def test_lower_central_series():
    G = SymmetricGroup.as_permutation_group(4)
    assert not G.is_nilpotent()