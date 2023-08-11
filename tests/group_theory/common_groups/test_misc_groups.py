"""Unit tests for misc groups (e.g. polyhedral groups, Klein 4 group, etc.)"""

import pytest
import quaternionic

from noetherpy.group_theory.group_elements import QuaternionElement

from noetherpy.group_theory.common_groups.misc_groups import (
    dicyclic_group,
    dihedral_group,
    klein_four_group,
    modular_maximal_cyclic_group,
    quasidihedral_group,
    quaternion_group,
)


def test_no_non_positive_dihedral_group_orders() -> None:
    with pytest.raises(ValueError):
        dihedral_group(0)
    with pytest.raises(ValueError):
        dihedral_group(-1)


def test_no_non_positive_dicyclic_group_orders() -> None:
    with pytest.raises(ValueError):
        dicyclic_group(0)
    with pytest.raises(ValueError):
        dicyclic_group(-1)


def test_no_non_positive_quasidihedral_group_orders() -> None:
    with pytest.raises(ValueError):
        quasidihedral_group(0)
    with pytest.raises(ValueError):
        quasidihedral_group(-1)


def test_no_non_positive_modular_maximal_cyclic_group_orders() -> None:
    with pytest.raises(ValueError):
        modular_maximal_cyclic_group(0)
    with pytest.raises(ValueError):
        modular_maximal_cyclic_group(-1)


klein_4 = klein_four_group()
klein_4_permutation = klein_four_group(representation="permutation")
klein_4_matrix = klein_four_group(representation="matrix")

quaternions = quaternion_group()

dihedral_14 = dihedral_group(7)
dihedral_14_permutation = dihedral_group(sides=7, representation="permutation")
dihedral_14_matrix = dihedral_group(sides=7, representation="matrix")

dihedral_16 = dihedral_group(8)
dihedral_16_permutation = dihedral_group(sides=8, representation="permutation")
dihedral_16_matrix = dihedral_group(sides=8, representation="matrix")


def test_klein_4_group_order_is_4() -> None:
    assert klein_4.order == 4
    assert klein_4_permutation.order == 4
    assert klein_4_matrix.order == 4


def test_klein_4_group_is_abelian() -> None:
    assert klein_4.is_abelian() is True
    assert klein_4_permutation.is_abelian() is True
    assert klein_4_matrix.is_abelian() is True


def test_klein_4_group_equals_own_center() -> None:
    assert klein_4.center == klein_4
    assert klein_4_permutation.center == klein_4_permutation
    assert klein_4_matrix.center == klein_4_matrix


def test_quaternion_group_order_is_8() -> None:
    assert quaternions.order == 8


def test_quaternion_group_is_non_abelian() -> None:
    assert quaternions.is_abelian() is False


def test_center_of_quaternion_group() -> None:
    assert quaternions.center == quaternions.subgroup_generated_by(
        [-QuaternionElement(quaternionic.array(1, 0, 0, 0))]
    )


def test_dihedral_group_order():
    assert dihedral_14.order == 14
    assert dihedral_14_permutation.order == 14
    assert dihedral_14_matrix.order == 14

    assert dihedral_16.order == 16
    assert dihedral_16_permutation.order == 16
    assert dihedral_16_matrix.order == 16


def test_dihedral_group_center():
    assert dihedral_14.center.is_trivial()
    assert dihedral_14_permutation.center.is_trivial()
    assert dihedral_14_matrix.center.is_trivial()

    assert dihedral_16.center.order == 2
    assert dihedral_16_permutation.center.order == 2
    assert dihedral_16_matrix.center.order == 2


def test_dicyclic_group_order() -> None:
    assert dicyclic_group(7).order == 4 * 7
    assert dicyclic_group(8).order == 4 * 8
    assert dicyclic_group(9).order == 4 * 9


def test_quasidihedral_group_order() -> None:
    assert quasidihedral_group(2).order == 2**2
    assert quasidihedral_group(3).order == 2**3
    assert quasidihedral_group(4).order == 2**4


def test_modular_maximal_cyclic_group_order() -> None:
    assert modular_maximal_cyclic_group(2).order == 2**2
    assert modular_maximal_cyclic_group(3).order == 2**3
    assert modular_maximal_cyclic_group(4).order == 2**4
