"""Tests for group homomorphisms."""


from noetherpy.group_theory.group_elements import CyclicGroupElement, GroupElement
from noetherpy.group_theory.homomorphisms import Homomorphism
from noetherpy.group_theory.common_groups.cyclic_groups import cyclic_group
from noetherpy.group_theory.common_groups.finite_abelian_groups import (
    finite_abelian_group_from_order_power_dict,
)
from noetherpy.group_theory.common_groups.misc_groups import dihedral_group
from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group

G = symmetric_group(4)
g = G[2]


def test_is_identity() -> None:
    def identity_map(x: GroupElement) -> GroupElement:
        return x

    id_map = Homomorphism(G, identity_map, G)

    assert id_map.is_identity()


def test_is_trivial() -> None:
    def identity_map(x: GroupElement) -> GroupElement:
        del x
        return G.identity

    trivial_map = Homomorphism(G, identity_map, G)

    assert trivial_map.is_trivial()


def test_identity_map_is_valid_automorphism() -> None:
    def identity_map(x: GroupElement) -> GroupElement:
        return x

    id_map = Homomorphism(G, identity_map, G)

    assert id_map.is_iso


def test_conjugation_is_valid_automorphism() -> None:
    def conjugation(x: GroupElement) -> GroupElement:
        return g * x * (~g)

    conj_auto = Homomorphism(G, conjugation, G)

    assert conj_auto.is_iso


def test_left_multiplication_is_not_valid_homomorphism() -> None:
    def left_multiplication(x: GroupElement) -> GroupElement:
        return g * x

    left_mult = Homomorphism(G, left_multiplication, G)

    assert ~left_mult.validate_homomorphism()


def test_squaring_is_valid_homomorphism_for_abelian_groups() -> None:
    A = finite_abelian_group_from_order_power_dict({3: 1, 6: 1, 7: 1})

    def squaring_map(x: GroupElement) -> GroupElement:
        return x * x

    assert ~Homomorphism(A, squaring_map, A).validate_homomorphism()


def test_inversion_is_valid_homomorphism_for_abelian_groups() -> None:
    A = finite_abelian_group_from_order_power_dict({3: 1, 6: 1, 7: 1})

    def squaring_map(x: GroupElement) -> GroupElement:
        return ~x

    assert ~Homomorphism(A, squaring_map, A).validate_homomorphism()


D = dihedral_group(8)
C = cyclic_group(8)


# this is the homomorphism which exhibits the subgroup
# of rotations as the kernel of a homomorphism
def rotations_to_kernel(x):
    return CyclicGroupElement(8, 4 * x.exponents[1])


rtk = Homomorphism(D, rotations_to_kernel, C)


def test_kernel_is_normal_subgroup() -> None:
    assert rtk.validate_homomorphism()
    assert rtk.kernel.is_normal


def test_image_is_subgroup() -> None:
    assert rtk.image == C.subgroup_generated_by([CyclicGroupElement(8, 4)])
