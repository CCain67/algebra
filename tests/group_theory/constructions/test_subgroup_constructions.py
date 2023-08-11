"""Unit tests for subgroup constructions"""

import pytest

from noetherpy.group_theory.groups import Subgroup
from noetherpy.group_theory.group_elements import DihedralGroupElement, Permutation
from noetherpy.group_theory.common_groups.permutation_groups import (
    alternating_group,
    symmetric_group,
)
from noetherpy.group_theory.common_groups.misc_groups import dihedral_group

from noetherpy.group_theory.constructions.subgroup_constructions import (
    centralizer,
    conjugate_subgroup,
    normal_closure,
    normal_core,
    normalizer,
)

sym_4 = symmetric_group(4)
alt_4 = Subgroup(alternating_group(4), sym_4)
sym_sub = sym_4.subgroup_generated_by([Permutation([1, 2, 3, 0])])

dih_7 = dihedral_group(7)
rot_7 = dih_7.subgroup_generated_by([DihedralGroupElement((1, 0), 7)])
flip_7 = dih_7.subgroup_generated_by([DihedralGroupElement((0, 1), 7)])


def test_centralizer_parent_group_value_error() -> None:
    with pytest.raises(ValueError):
        centralizer(flip_7, sym_4)


def test_normalizer_parent_group_value_error() -> None:
    with pytest.raises(ValueError):
        normalizer(flip_7, sym_4)


def test_centralizer_contained_in_normalizer() -> None:
    assert centralizer(sym_sub, sym_4) < normalizer(sym_sub, sym_4)
    assert centralizer(flip_7, dih_7) < normalizer(flip_7, dih_7)


def test_centralizer_is_normal_subgroup_of_normalizer() -> None:
    assert Subgroup(
        centralizer(sym_sub, sym_4).elements, normalizer(sym_sub, sym_4)
    ).is_normal
    assert Subgroup(
        centralizer(flip_7, dih_7).elements, normalizer(flip_7, dih_7)
    ).is_normal


def test_normalizer_of_normal_subgroup_is_entire_group() -> None:
    assert normalizer(alt_4, sym_4) == sym_4
    assert normalizer(rot_7, dih_7) == dih_7


def test_subgroup_is_normal_inside_of_normalizer() -> None:
    assert Subgroup(sym_sub.elements, normalizer(sym_sub, sym_4)).is_normal
    assert Subgroup(flip_7.elements, normalizer(flip_7, dih_7)).is_normal
