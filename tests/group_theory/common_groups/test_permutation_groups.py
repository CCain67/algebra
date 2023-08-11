"""Unit tests for permutation groups"""
import pytest

from noetherpy.group_theory.groups import Subgroup
from noetherpy.group_theory.common_groups.permutation_groups import (
    alternating_group,
    symmetric_group,
)


def test_no_non_positive_symmetric_group_orders() -> None:
    with pytest.raises(ValueError):
        symmetric_group(0)
    with pytest.raises(ValueError):
        symmetric_group(-1)


def test_no_non_positive_alternating_group_orders() -> None:
    with pytest.raises(ValueError):
        alternating_group(0)
    with pytest.raises(ValueError):
        alternating_group(-1)


sym_2 = symmetric_group(2)
sym_3 = symmetric_group(3)
sym_4 = symmetric_group(4)

sym_2_matrix = symmetric_group(N=2, representation="matrix")
sym_3_matrix = symmetric_group(N=3, representation="matrix")
sym_4_matrix = symmetric_group(N=4, representation="matrix")

alt_2 = Subgroup(alternating_group(2).elements, sym_2)
alt_3 = Subgroup(alternating_group(3).elements, sym_3)
alt_4 = Subgroup(alternating_group(4).elements, sym_4)

alt_2_matrix = Subgroup(
    alternating_group(N=2, representation="matrix").elements, sym_2_matrix
)
alt_3_matrix = Subgroup(
    alternating_group(N=3, representation="matrix").elements, sym_3_matrix
)
alt_4_matrix = Subgroup(
    alternating_group(N=4, representation="matrix").elements, sym_4_matrix
)


def test_symmetric_group_order():
    assert sym_2.order == 2
    assert sym_3.order == 6
    assert sym_4.order == 24

    assert sym_2_matrix.order == 2
    assert sym_3_matrix.order == 6
    assert sym_4_matrix.order == 24


def test_alternating_group_order():
    assert alt_2.order == 1
    assert alt_3.order == 3
    assert alt_4.order == 12

    assert alt_2_matrix.order == 1
    assert alt_3_matrix.order == 3
    assert alt_4_matrix.order == 12


def test_symmetric_group_abelian_cases():
    assert sym_2.is_abelian() is True
    assert sym_3.is_abelian() is False
    assert sym_4.is_abelian() is False

    assert sym_2_matrix.is_abelian() is True
    assert sym_3_matrix.is_abelian() is False
    assert sym_4_matrix.is_abelian() is False


def test_alternating_group_abelian_cases():
    assert alt_2.is_abelian() is True
    assert alt_3.is_abelian() is True
    assert alt_4.is_abelian() is False

    assert alt_2_matrix.is_abelian() is True
    assert alt_3_matrix.is_abelian() is True
    assert alt_4_matrix.is_abelian() is False


def test_symmetric_group_trivial_center_cases():
    assert sym_2.center.is_trivial() is False
    assert sym_3.center.is_trivial() is True
    assert sym_4.center.is_trivial() is True

    assert sym_2_matrix.center.is_trivial() is False
    assert sym_3_matrix.center.is_trivial() is True
    assert sym_4_matrix.center.is_trivial() is True


def test_alternating_group_trivial_center_cases():
    assert alt_2.center.is_trivial() is True
    assert alt_3.center.is_trivial() is False
    assert alt_4.center.is_trivial() is True

    assert alt_2_matrix.center.is_trivial() is True
    assert alt_3_matrix.center.is_trivial() is False
    assert alt_4_matrix.center.is_trivial() is True


def test_symmetric_group_conjugacy_classes() -> None:
    assert len(sym_2.conjugacy_classes) == 2
    assert len(sym_3.conjugacy_classes) == 3
    assert len(sym_4.conjugacy_classes) == 5
    assert len(sym_2_matrix.conjugacy_classes) == 2
    assert len(sym_3_matrix.conjugacy_classes) == 3
    assert len(sym_4_matrix.conjugacy_classes) == 5


def test_alternating_group_conjugacy_classes() -> None:
    assert len(alt_2.conjugacy_classes) == 1
    assert len(alt_3.conjugacy_classes) == 3
    assert len(alt_4.conjugacy_classes) == 4
    assert len(alt_2_matrix.conjugacy_classes) == 1
    assert len(alt_3_matrix.conjugacy_classes) == 3
    assert len(alt_4_matrix.conjugacy_classes) == 4


def test_alternating_group_inclusion_into_symmetric_group() -> None:
    assert alt_2.validate_inclusion()
    assert alt_2_matrix.validate_inclusion()
    assert alt_3.validate_inclusion()
    assert alt_3_matrix.validate_inclusion()
    assert alt_4.validate_inclusion()
    assert alt_4_matrix.validate_inclusion()


def test_alternating_group_is_subgroup_of_symmetric_group() -> None:
    assert alt_2.validate_subgroup()
    assert alt_2_matrix.validate_subgroup()
    assert alt_3.validate_subgroup()
    assert alt_3_matrix.validate_subgroup()
    assert alt_4.validate_subgroup()
    assert alt_4_matrix.validate_subgroup()


def test_alternating_group_is_normal_subgroup_of_symmetric_group() -> None:
    assert alt_2.is_normal
    assert alt_2_matrix.is_normal
    assert alt_3.is_normal
    assert alt_3_matrix.is_normal
    assert alt_4.is_normal
    assert alt_4_matrix.is_normal
