import itertools

from group_theory.group_elements import (
    Permutation,
)
from group_theory.groups import Group


def symmetric_group(N: int, repr: str = "permutation") -> Group:
    if N < 1:
        raise ValueError(
            "the symmetricing group on N letters requires N to be at least 1"
        )
    if repr not in ["permutation", "matrix"]:
        raise ValueError('repr must be one of: "permutation" or "matrix"')

    P = list(itertools.permutations(range(1, N + 1)))
    S_N = Group([Permutation({i + 1: p[i] for i in range(N)}) for p in P])

    t = {1: 2, 2: 1}
    for i in range(3, N + 1):
        t[i] = i
    c = {i: i + 1 for i in range(1, N)}
    c[N] = 1
    S_N.canonical_generators = list({Permutation(c), Permutation(t)})

    if repr == "matrix":
        S_N.elements = [P.to_matrix() for P in S_N.elements]
        S_N.canonical_generators = list(
            {Permutation(c).to_matrix(), Permutation(t).to_matrix()}
        )
    return S_N


def alternating_group(N: int, repr: str = "permutation") -> Group:
    if N < 1:
        raise ValueError(
            "the alternating group on N letters requires N to be at least 1"
        )
    if repr not in ["permutation", "matrix"]:
        raise ValueError('repr must be one of: "permutation" or "matrix"')

    P = list(itertools.permutations(range(1, N + 1)))
    S_N = [Permutation({i + 1: p[i] for i in range(N)}) for p in P]
    A_N = Group([P for P in S_N if P.sign == 1])

    if N >= 3:
        A_N.canonical_generators = [
            Permutation({1: 2, 2: K, K: 1, **{i: i for i in range(3, N + 1) if i != K}})
            for K in range(3, N + 1)
        ]

    if repr == "matrix":
        A_N.elements = [P.to_matrix() for P in S_N.elements]
        A_N.canonical_generators = [P.to_matrix() for P in A_N.canonical_generators]
    return A_N
