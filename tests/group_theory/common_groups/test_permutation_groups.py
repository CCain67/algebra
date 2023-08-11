"""Unit tests for permutation groups"""
from noetherpy.group_theory.common_groups.permutation_groups import (
    alternating_group,
    symmetric_group,
)

sym_2 = symmetric_group(2)
sym_3 = symmetric_group(3)
sym_4 = symmetric_group(4)

sym_2_matrix = symmetric_group(N=2, representation="matrix")
sym_3_matrix = symmetric_group(N=3, representation="matrix")
sym_4_matrix = symmetric_group(N=4, representation="matrix")

alt_2 = alternating_group(2)
alt_3 = alternating_group(3)
alt_4 = alternating_group(4)

alt_2_matrix = alternating_group(N=2, representation="matrix")
alt_3_matrix = alternating_group(N=3, representation="matrix")
alt_4_matrix = alternating_group(N=4, representation="matrix")


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
