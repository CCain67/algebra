"""This module defines the functions which construct the symmetric and alternating groups."""

import itertools

from ..group_elements import (
    Permutation,
)
from ..groups import Group


def symmetric_group(N: int, representation: str = "permutation") -> Group:
    """Constructs the symmetric group on N symbols.

    Args:
        N (int): number of symbols
        representation (str, optional): Representation of the group elements.
        Defaults to "permutation".

    Raises:
        ValueError: Raised if N < 1.
        ValueError: Raised if a bad representation option is passed.

    Returns:
        Group: the symmetric group on N syumbols.
    """
    if N < 1:
        raise ValueError(
            "the symmetricing group on N letters requires N to be at least 1"
        )
    if representation not in ["permutation", "matrix"]:
        raise ValueError('repr must be one of: "permutation" or "matrix"')

    all_permutations = list(itertools.permutations(range(1, N + 1)))
    sym_group = Group(
        [Permutation({i + 1: p[i] for i in range(N)}) for p in all_permutations]
    )

    transposition = {1: 2, 2: 1}
    for i in range(3, N + 1):
        transposition[i] = i
    cycle = {i: i + 1 for i in range(1, N)}
    cycle[N] = 1
    sym_group.canonical_generators = list(
        {Permutation(cycle), Permutation(transposition)}
    )

    if representation == "matrix":
        sym_group.elements = [P.to_matrix() for P in sym_group.elements]
        sym_group.canonical_generators = list(
            {Permutation(cycle).to_matrix(), Permutation(transposition).to_matrix()}
        )
    return sym_group


def alternating_group(N: int, representation: str = "permutation") -> Group:
    """Constructs the alternating group on N symbols.

    Args:
        N (int): number of symbols
        representation (str, optional): Representation of the group elements.
        Defaults to "permutation".

    Raises:
        ValueError: Raised if N < 1.
        ValueError: Raised if a bad representation option is passed.

    Returns:
        Group: the alternating group on N syumbols.
    """
    if N < 1:
        raise ValueError(
            "the alternating group on N letters requires N to be at least 1"
        )
    if representation not in ["permutation", "matrix"]:
        raise ValueError('repr must be one of: "permutation" or "matrix"')

    if N == 2:
        canonical_generators = [Permutation({1: 2, 2: 1})]
        identity = Permutation({1: 1, 2: 2})

    if N >= 3:
        canonical_generators = [
            Permutation({1: 2, 2: K, K: 1, **{i: i for i in range(3, N + 1) if i != K}})
            for K in range(3, N + 1)
        ]
        identity = Permutation({i: i for i in range(1, N + 1)})

    if representation == "matrix":
        canonical_generators = [P.to_matrix() for P in canonical_generators]
        identity = identity.to_matrix()
    return Group.from_generators(canonical_generators, identity)
