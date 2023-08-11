"""This module defines functions for constructing the cyclic groups, 
the Klein 4 group, and dihedral groups."""

import quaternionic

from ..group_elements import (
    DicyclicGroupElement,
    DihedralGroupElement,
    ModularMaximalCyclicGroupElement,
    Permutation,
    QuasidihedralGroupElement,
    QuaternionElement,
)
from ..groups import Group
from ..common_groups.cyclic_groups import cyclic_group


def klein_four_group(representation: str = "residue") -> Group:
    """Constructs the Klein four group.

    Args:
        repr (str, optional): Representation of the group. Defaults to "residue".

    Raises:
        ValueError: This error is raised if a bad representation option is passed.

    Returns:
        Group: The Klein four group.
    """
    if representation not in ["residue", "permutation", "matrix"]:
        raise ValueError('repr must be one of: "residue", "permutation" or "matrix"')

    if representation == "residue":
        cyclic_order_2 = cyclic_group(2)
        return cyclic_order_2 * cyclic_order_2
    if representation == "permutation":
        return Group.from_generators(
            [
                Permutation([1, 0, 2, 3]),
                Permutation([0, 1, 3, 2]),
            ],
            Permutation(list(range(4))),
        )
    return Group.from_generators(
        [
            Permutation([1, 0, 2, 3]).to_matrix(),
            Permutation([0, 1, 3, 2]).to_matrix(),
        ],
        Permutation(list(range(4))).to_matrix(),
    )


def quaternion_group(representation: str = "quaternion") -> Group:
    """Constructs the quaternion group.

    Args:
        representation (str, optional): choice of representation for
        the quaternion elements. Defaults to "quaternion".

    Raises:
        ValueError: this is raised if a bad representation option is passed

    Returns:
        Group: the quaternion group of order 8.
    """
    if representation not in ["quaternion", "permutation", "matrix"]:
        raise ValueError('repr must be one of: "quaternion", "permutation" or "matrix"')

    if representation != "quaternion":
        return NotImplemented

    one = QuaternionElement(quaternionic.array(1, 0, 0, 0))
    i = QuaternionElement(quaternionic.i)
    j = QuaternionElement(quaternionic.j)
    k = QuaternionElement(quaternionic.k)

    return Group([one, i, j, k, -one, -i, -j, -k])


def dihedral_group(sides: int, representation: str = "dihedral") -> Group:
    """Constructs the dihedral group - the symmetry group of the n-gon.

    Args:
        sides (int): the number of sides of the n-gon to form
        the symmetry group of.
        representation (str, optional): Representation of the group. Defaults to "permutation".

    Raises:
        ValueError: This error is raised if a bad representation option is passed.

    Returns:
        Group: The dihedral group of order 2*sides.
    """
    if representation not in ["dihedral", "permutation", "matrix"]:
        raise ValueError('repr must be one of: "permutation" or "matrix"')

    r = DihedralGroupElement((1, 0), sides)
    s = DihedralGroupElement((0, 1), sides)

    if representation == "permutation":
        r = r.to_permutation()
        s = s.to_permutation()
    elif representation == "matrix":
        r = r.to_matrix()
        s = s.to_matrix()
    dih_n = Group([(r**k) * (s**i) for k in range(sides) for i in range(2)])
    dih_n.canonical_generators = [r, s]
    return dih_n


def dicyclic_group(n: int) -> Group:
    """Constructs the nth dicyclic group.

    Args:
        n (int): The parameter defining the dicyclic group.

    Returns:
        Group: the nth dicyclic group.
    """
    a = DicyclicGroupElement((1, 0), n)
    x = DicyclicGroupElement((0, 1), n)
    dic_n = Group(
        [DicyclicGroupElement((k, i), n) for k in range(2 * n) for i in range(2)]
    )
    dic_n.canonical_generators = [a, x]
    return dic_n


def quasidihedral_group(n: int) -> Group:
    """Constructs the nth quasidihedral group.

    Args:
        n (int): The parameter defining the quasidihedral group.

    Returns:
        Group: the nth quasidihedral group.
    """
    r = QuasidihedralGroupElement((1, 0), n)
    s = QuasidihedralGroupElement((0, 1), n)
    qdih_n = Group(
        [
            QuasidihedralGroupElement((k, i), n)
            for k in range(2 ** (n - 1))
            for i in range(2)
        ]
    )
    qdih_n.canonical_generators = [r, s]
    return qdih_n


def modular_maximal_cyclic_group(n: int) -> Group:
    """Constructs the nth modular maximal-cyclic group.

    Args:
        n (int): The parameter defining the modular maximal-cyclic group.

    Returns:
        Group: the nth modular maximal-cyclic group.
    """
    r = ModularMaximalCyclicGroupElement((1, 0), n)
    s = ModularMaximalCyclicGroupElement((0, 1), n)
    mmc_n = Group(
        [
            ModularMaximalCyclicGroupElement((k, i), n)
            for k in range(2 ** (n - 1))
            for i in range(2)
        ]
    )
    mmc_n.canonical_generators = [r, s]
    return mmc_n
