"""Unit tests for cyclic groups"""

import pytest

from noetherpy.group_theory.common_groups.cyclic_groups import (
    cyclic_group,
)

cyclic_group_1 = cyclic_group(1)
cyclic_group_6 = cyclic_group(6)
cyclic_group_7 = cyclic_group(7)

cyclic_group_1_permutation = cyclic_group(N=1, representation="permutation")
cyclic_group_6_permutation = cyclic_group(N=6, representation="permutation")
cyclic_group_7_permutation = cyclic_group(N=7, representation="permutation")

cyclic_group_1_matrix = cyclic_group(N=1, representation="matrix")
cyclic_group_6_matrix = cyclic_group(N=6, representation="matrix")
cyclic_group_7_matrix = cyclic_group(N=7, representation="matrix")


def test_no_non_positive_cyclic_group_orders() -> None:
    with pytest.raises(ValueError):
        cyclic_group(0)
    with pytest.raises(ValueError):
        cyclic_group(-1)


def test_cyclic_group_triviality():
    assert cyclic_group_1.is_trivial() is True
    assert cyclic_group_6.is_trivial() is False
    assert cyclic_group_7.is_trivial() is False

    assert cyclic_group_1_permutation.is_trivial() is True
    assert cyclic_group_6_permutation.is_trivial() is False
    assert cyclic_group_7_permutation.is_trivial() is False

    assert cyclic_group_1_matrix.is_trivial() is True
    assert cyclic_group_6_matrix.is_trivial() is False
    assert cyclic_group_7_matrix.is_trivial() is False


def test_cyclic_group_abelian():
    assert cyclic_group_1.is_abelian() is True
    assert cyclic_group_6.is_abelian() is True
    assert cyclic_group_7.is_abelian() is True

    assert cyclic_group_1_permutation.is_abelian() is True
    assert cyclic_group_6_permutation.is_abelian() is True
    assert cyclic_group_7_permutation.is_abelian() is True

    assert cyclic_group_1_matrix.is_abelian() is True
    assert cyclic_group_6_matrix.is_abelian() is True
    assert cyclic_group_7_matrix.is_abelian() is True


def test_cyclic_group_order():
    assert cyclic_group_1.order == 1
    assert cyclic_group_6.order == 6
    assert cyclic_group_7.order == 7

    assert cyclic_group_1_permutation.order == 1
    assert cyclic_group_6_permutation.order == 6
    assert cyclic_group_7_permutation.order == 7

    assert cyclic_group_1_matrix.order == 1
    assert cyclic_group_6_matrix.order == 6
    assert cyclic_group_7_matrix.order == 7


def test_cyclic_group_conjugacy_classes():
    assert len(cyclic_group_1.conjugacy_classes) == 1
    assert len(cyclic_group_6.conjugacy_classes) == 6
    assert len(cyclic_group_7.conjugacy_classes) == 7

    assert len(cyclic_group_1_permutation.conjugacy_classes) == 1
    assert len(cyclic_group_6_permutation.conjugacy_classes) == 6
    assert len(cyclic_group_7_permutation.conjugacy_classes) == 7

    assert len(cyclic_group_1_matrix.conjugacy_classes) == 1
    assert len(cyclic_group_6_matrix.conjugacy_classes) == 6
    assert len(cyclic_group_7_matrix.conjugacy_classes) == 7


def test_cyclic_group_equals_own_center():
    assert cyclic_group_1.center == cyclic_group_1
    assert cyclic_group_6.center == cyclic_group_6
    assert cyclic_group_7.center == cyclic_group_7

    assert cyclic_group_1_permutation.center == cyclic_group_1_permutation
    assert cyclic_group_6_permutation.center == cyclic_group_6_permutation
    assert cyclic_group_7_permutation.center == cyclic_group_7_permutation

    assert cyclic_group_1_matrix.center == cyclic_group_1_matrix
    assert cyclic_group_6_matrix.center == cyclic_group_6_matrix
    assert cyclic_group_7_matrix.center == cyclic_group_7_matrix
