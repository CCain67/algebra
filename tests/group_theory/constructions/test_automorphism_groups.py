"""Unit tests for automorphism group computations"""

from noetherpy.group_theory.constructions.automorphism_groups import Aut, Inn, Out

from noetherpy.group_theory.common_groups.permutation_groups import (
    alternating_group,
    symmetric_group,
)
from noetherpy.group_theory.common_groups.misc_groups import (
    dihedral_group,
    klein_four_group,
)


dih_7 = dihedral_group(7)
dih_8 = dihedral_group(8)

sym_5 = symmetric_group(5)
alt_4 = alternating_group(4)
klein_4 = klein_four_group()


def test_inner_automorphism_group_has_correct_order() -> None:
    assert Inn(klein_4).order == 1
    assert Inn(alt_4).order == 12
    assert Inn(sym_5).order == 120


def test_automorphism_group_has_correct_order() -> None:
    assert Aut(klein_4).order == 6
    assert Aut(alt_4).order == 24
    assert Aut(sym_5).order == 120


def test_outer_automorphism_group_has_correct_order() -> None:
    assert Out(klein_4).order == 6
    assert Out(alt_4).order == 2
    assert Out(sym_5).order == 1
