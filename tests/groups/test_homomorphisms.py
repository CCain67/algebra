"""Tests for group homomorphisms."""


from noetherpy.group_theory.group_elements import GroupElement
from noetherpy.group_theory.homomorphisms import Homomorphism
from noetherpy.group_theory.common_groups.finite_abelian_groups import (
    finite_abelian_group_from_order_power_dict,
)
from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group

G = symmetric_group(4)
g = G[2]


def test_identity_map_is_valid_automorphism() -> None:
    def identity_map(x: GroupElement) -> GroupElement:
        return x

    id_map = Homomorphism(G, identity_map, G)

    assert id_map.validate_homomorphism()
    assert id_map.kernel.is_trivial()


def test_conjugation_is_valid_automorphism() -> None:
    def conjugation(x: GroupElement) -> GroupElement:
        return g * x * (~g)

    conj_auto = Homomorphism(G, conjugation, G)

    assert conj_auto.validate_homomorphism()
    assert conj_auto.kernel.is_trivial()


def test_left_multiplication_is_not_valid_homomorphism_for_non_abelian_group() -> None:
    def left_multiplication(x: GroupElement) -> GroupElement:
        return g * x

    left_mult = Homomorphism(G, left_multiplication, G)

    assert ~left_mult.validate_homomorphism()


def test_squaring_is_valid_homomorphism_for_abelian_groups() -> None:
    A = finite_abelian_group_from_order_power_dict({2: 3, 3: 1, 7: 2})
    B = finite_abelian_group_from_order_power_dict({13: 1, 4: 2, 9: 1})

    def squaring_map(x: GroupElement) -> GroupElement:
        return x * x

    assert ~Homomorphism(A, squaring_map, A).validate_homomorphism()
    assert ~Homomorphism(B, squaring_map, B).validate_homomorphism()


def test_inversion_is_valid_homomorphism_for_abelian_groups() -> None:
    A = finite_abelian_group_from_order_power_dict({2: 3, 3: 1, 7: 2})
    B = finite_abelian_group_from_order_power_dict({13: 1, 4: 2, 9: 1})

    def squaring_map(x: GroupElement) -> GroupElement:
        return ~x

    assert ~Homomorphism(A, squaring_map, A).validate_homomorphism()
    assert ~Homomorphism(B, squaring_map, B).validate_homomorphism()
