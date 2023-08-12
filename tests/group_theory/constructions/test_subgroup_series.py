"""Unit tests for subgroup series"""

from noetherpy.group_theory.common_groups.permutation_groups import (
    alternating_group,
    symmetric_group,
)
from noetherpy.group_theory.common_groups.misc_groups import (
    dihedral_group,
    quaternion_group,
)
from noetherpy.group_theory.constructions.subgroup_series import (
    derived_series,
    is_solvable,
    is_nilpotent,
    lower_central_series,
    upper_central_series,
)


sym_4 = symmetric_group(4)
sym_5 = symmetric_group(5)
alt_5 = alternating_group(5)

dih_5 = dihedral_group(5)
dih_6 = dihedral_group(6)
quat = quaternion_group()


def test_derived_series_has_correct_length() -> None:
    assert len(derived_series(sym_4)) == 4
    assert len(derived_series(dih_5)) == 3
    assert len(derived_series(dih_6)) == 3
    assert len(derived_series(quat)) == 3


def test_lower_central_series_has_correct_length() -> None:
    assert len(lower_central_series(sym_4)) == 2
    assert len(lower_central_series(dih_5)) == 2
    assert len(lower_central_series(dih_6)) == 2
    assert len(lower_central_series(quat)) == 3


def test_upper_central_series_has_correct_length() -> None:
    assert len(upper_central_series(sym_4)) == 1
    assert len(upper_central_series(dih_5)) == 1
    assert len(upper_central_series(dih_6)) == 2
    assert len(upper_central_series(quat)) == 3


def test_is_solvable() -> None:
    assert is_solvable(quat)
    assert is_solvable(sym_4)
    assert not is_solvable(sym_5)
    assert not is_solvable(alt_5)


def test_is_nilpotent() -> None:
    assert is_nilpotent(quat)
    assert not is_nilpotent(sym_4)
