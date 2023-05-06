"""This module defines functions for constructing the cyclic groups, 
the Klein 4 group, and dihedral groups."""

import numpy

from group_theory.group_elements import (
    AdditiveResidueClass,
    Matrix,
    Permutation,
)
from group_theory.groups import Group


def cyclic_group(N: int, representation: str = "residue") -> Group:
    """Constructs the cyclic group of order N.

    Args:
        N (int): order of the cyclic group to be constructed
        representation (str, optional): Representation of the cyclic group elements.
        Defaults to "residue".

    Raises:
        ValueError: Raised if a bad representation option is passed.

    Returns:
        Group: the cyclic group of order N.
    """
    if representation not in ["residue", "permutation"]:
        raise ValueError('repr must be one of: "residue" or "permutation"')

    if representation == "residue":
        g = AdditiveResidueClass(1, N)
    elif representation == "permutation":
        generator = {i: i + 1 for i in range(1, N)}
        generator[N] = 1
        g = Permutation(generator)
    elif representation == "matrix":
        generator = {i: i + 1 for i in range(1, N)}
        generator[N] = 1
        matrix_generator = numpy.zeros((N, N))
        for k in generator.keys():
            matrix_generator[k - 1, generator[k] - 1] = 1
        g = Matrix(matrix_generator, 2, 1)

    cyc_group = Group([g**j for j in range(N)])
    cyc_group.canonical_generators = [g]
    return cyc_group


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


def dihedral_group(sides: int, representation: str = "permutation") -> Group:
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
    if representation not in ["permutation", "matrix"]:
        raise ValueError('repr must be one of: "permutation" or "matrix"')

    cycle = {**{i: i + 1 for i in range(1, sides)}, sides: 1}
    flip = {1: 1, **{2 + i: sides - i for i in range(sides - 1)}}

    if representation == "permutation":
        r = Permutation(cycle)
        s = Permutation(flip)
    elif representation == "matrix":
        r = Permutation(cycle).to_matrix()
        s = Permutation(flip).to_matrix()

    dihedral_group_elements = [
        (r**i) * (s**j) for i in range(sides) for j in range(2)
    ]
    dh_group = Group(dihedral_group_elements)
    dh_group.canonical_generators = [r, s]
    return dh_group
