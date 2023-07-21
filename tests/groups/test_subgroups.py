import pytest

from group_theory.groups import Subgroup

from group_theory.common_groups.permutation_groups import symmetric_group
from group_theory.common_groups.misc_groups import dihedral_group

from group_theory.constructions.subgroup_constructions import centralizer, normalizer
from group_theory.constructions.subgroup_series import is_solvable, is_nilpotent


def test_center_dihedral():
    assert dihedral_group(7).center.order == 1
    assert dihedral_group(8).center.order == 2


def test_normalizer():
    G = symmetric_group(5)
    S = G.subgroup_generated_by([G.canonical_generators[0]])
    N = normalizer(S, G)
    C = centralizer(S, G)
    S = Subgroup(S.elements, N)
    assert S.is_normal
    assert C < N


def test_derived_series():
    S_4 = symmetric_group(4)
    S_5 = symmetric_group(5)
    assert is_solvable(S_4)
    assert not is_solvable(S_5)


def test_lower_central_series():
    G = symmetric_group(4)
    assert not is_nilpotent(G)
