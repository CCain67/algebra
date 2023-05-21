"""This module defines functions for constructing the cyclic groups, 
the Klein 4 group, and dihedral groups."""

from group_theory.group_elements import (
    DicyclicGroupElement,
    DihedralGroupElement,
    QuasidihedralGroupElement,
)
from group_theory.groups import Group
from group_theory.common_groups.cyclic_groups import cyclic_group


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
    return dihedral_group(2, representation)


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
    a = QuasidihedralGroupElement((1, 0), n)
    x = QuasidihedralGroupElement((0, 1), n)
    qdih_n = Group(
        [
            QuasidihedralGroupElement((k, i), n)
            for k in range(2 ** (n - 1))
            for i in range(2)
        ]
    )
    qdih_n.canonical_generators = [a, x]
    return qdih_n
