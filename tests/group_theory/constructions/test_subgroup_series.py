"""Unit tests for subgroup series"""
from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group


from noetherpy.group_theory.constructions.subgroup_series import (
    is_solvable,
    is_nilpotent,
)


def test_derived_series():
    S_4 = symmetric_group(4)
    S_5 = symmetric_group(5)
    assert is_solvable(S_4)
    assert not is_solvable(S_5)


def test_lower_central_series():
    G = symmetric_group(4)
    assert not is_nilpotent(G)
