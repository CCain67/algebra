"""Unit tests for permutation groups"""
from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group


def test_symmetric_group_order():
    assert symmetric_group(4).order == 24
    assert symmetric_group(4, representation="matrix").order == 24


def test_symmetric_quotient_alternating():
    S = symmetric_group(4)
    A = S.commutator_subgroup
    assert (S / A).order == 2
