"""Unit tests for various types of products of groups"""

from noetherpy.group_theory.group_elements import DihedralGroupElement, GroupElement
from noetherpy.group_theory.homomorphisms import Homomorphism
from noetherpy.group_theory.constructions.products import (
    semidirect_product,
    SemidirectProductElement,
)

from noetherpy.group_theory.common_groups.cyclic_groups import cyclic_group
from noetherpy.group_theory.common_groups.misc_groups import dihedral_group

cyc_2 = cyclic_group(2)
cyc_7 = cyclic_group(7)

dih_14 = dihedral_group(7)


def inversion_map(g: GroupElement, x: GroupElement) -> GroupElement:
    assert g in cyc_2
    power_map = {0: 1, 1: -1}
    return x ** power_map[g.power]


# this is isomorphic to the dihedral group D_14
cyc_7_sdp_cyc_2 = semidirect_product(cyc_7, inversion_map, cyc_2)

# the isomorphism
isomorphism = Homomorphism.from_dict(
    cyc_7_sdp_cyc_2,
    {
        SemidirectProductElement(
            (cyc_7[i], cyc_2[j]), inversion_map
        ): DihedralGroupElement((i, j), 7)
        for i in range(7)
        for j in range(2)
    },
    dih_14,
)


def test_semidirect_product_has_correct_identity_element() -> None:
    assert cyc_7_sdp_cyc_2.identity == SemidirectProductElement(
        (cyc_7.identity, cyc_2.identity), inversion_map
    )


def test_semi_direct_product_has_correct_order() -> None:
    assert cyc_7_sdp_cyc_2.order == dih_14.order


def test_semidirect_product_gives_correct_group() -> None:
    assert isomorphism.check_iso()
