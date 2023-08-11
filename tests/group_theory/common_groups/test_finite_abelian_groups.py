"""Unit tests for finite abelian groups"""

from noetherpy.group_theory.common_groups.cyclic_groups import cyclic_group
from noetherpy.group_theory.common_groups.finite_abelian_groups import (
    finite_abelian_group_from_order_power_dict,
)


def test_finite_abelian_group_from_single_entry_order_power_dict_reduces_to_cyclic_group() -> None:
    assert finite_abelian_group_from_order_power_dict({7: 1}) == cyclic_group(7)
    assert finite_abelian_group_from_order_power_dict({12: 1}) == cyclic_group(12)


finite_abelian_group_12 = finite_abelian_group_from_order_power_dict({2: 2, 3: 1})
finite_abelian_group_30 = finite_abelian_group_from_order_power_dict({2: 1, 3: 1, 5: 1})

finite_abelian_group_12_permutation = finite_abelian_group_from_order_power_dict(
    order_power_dict={2: 2, 3: 1}, representation="permutation"
)
finite_abelian_group_30_permutation = finite_abelian_group_from_order_power_dict(
    order_power_dict={2: 1, 3: 1, 5: 1}, representation="permutation"
)

finite_abelian_group_12_matrix = finite_abelian_group_from_order_power_dict(
    order_power_dict={2: 2, 3: 1}, representation="matrix"
)
finite_abelian_group_30_matrix = finite_abelian_group_from_order_power_dict(
    order_power_dict={2: 1, 3: 1, 5: 1}, representation="matrix"
)


def test_finite_abelian_group_order() -> None:
    assert finite_abelian_group_12.order == 12
    assert finite_abelian_group_12_permutation.order == 12
    assert finite_abelian_group_12_matrix.order == 12

    assert finite_abelian_group_30.order == 30
    assert finite_abelian_group_30_permutation.order == 30
    assert finite_abelian_group_30_matrix.order == 30


def test_finite_abelian_group_is_abelian() -> None:
    assert finite_abelian_group_12.is_abelian() is True
    assert finite_abelian_group_12_permutation.is_abelian() is True
    assert finite_abelian_group_12_matrix.is_abelian() is True

    assert finite_abelian_group_30.is_abelian() is True
    assert finite_abelian_group_30_permutation.is_abelian() is True
    assert finite_abelian_group_30_matrix.is_abelian() is True

def test_finite_abelian_group_conjugacy_classes() -> None:
    assert len(finite_abelian_group_12.conjugacy_classes) == 12
    assert len(finite_abelian_group_12_permutation.conjugacy_classes) == 12
    assert len(finite_abelian_group_12_matrix.conjugacy_classes) == 12

    assert len(finite_abelian_group_30.conjugacy_classes) == 30
    assert len(finite_abelian_group_30_permutation.conjugacy_classes) == 30
    assert len(finite_abelian_group_30_matrix.conjugacy_classes) == 30

def test_finite_abelian_group_equals_own_center() -> None:
    assert finite_abelian_group_12.center == finite_abelian_group_12
    assert finite_abelian_group_12_permutation.center == finite_abelian_group_12_permutation
    assert finite_abelian_group_12_matrix.center == finite_abelian_group_12_matrix

    assert finite_abelian_group_30.center == finite_abelian_group_30
    assert finite_abelian_group_30_permutation.center == finite_abelian_group_30_permutation
    assert finite_abelian_group_30_matrix.center == finite_abelian_group_30_matrix
