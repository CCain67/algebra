"""Tests to guarantee quotient groups are computed correctly."""
from group_theory.common_groups.permutation_groups import symmetric_group


def test_symmetric_quotient_alternating():
    S = symmetric_group(4)
    A = S.commutator_subgroup
    assert (S / A).order == 2
