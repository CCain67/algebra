"""Tests to check if constructed groups contain the correct number of elements."""
from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group
from noetherpy.group_theory.common_groups.matrix_groups import (
    GL,
    heisenberg_group,
    O,
    SL,
)
from noetherpy.group_theory.common_groups.misc_groups import (
    cyclic_group,
    dihedral_group,
)


def test_cyclic_group_order():
    assert cyclic_group(31).order == 31


def test_symmetric_group_order():
    assert symmetric_group(4).order == 24
    assert symmetric_group(4, representation="matrix").order == 24


def test_dihedral_group():
    assert dihedral_group(7).order == 14


def get_GLnq_order(n, q):
    k = 0
    prod = 1
    while k <= n - 1:
        prod *= q**n - q**k
        k += 1
    return prod


def get_SLnq_order(n, q):
    return get_GLnq_order(n, q) / (q - 1)


def negative_one_is_square(q):
    if q % 2 == 0:
        # since -1=1 in a field of char=2, -1 is trivially a square
        return True
    return q % 4 == 1


# currently only for characteristic != 2
# see here: https://math.stackexchange.com/questions/564670/order-of-orthogonal-groups-over-finite-field
def get_Onq_order(n, q):
    if q % 2 == 0:
        return NotImplemented
    if n % 2 == 1:
        k = (n - 1) / 2
        i = 0
        prod = 1
        while i <= k - 1:
            prod *= q ** (2 * k) - q ** (2 * i)
            i += 1
        return 2 * (q**k) * prod
    k = n / 2
    i = 1
    prod = 1
    while i < k:
        prod *= q ** (2 * k) - q ** (2 * i)
        i += 1
    if negative_one_is_square(q):
        return 2 * (q**k - 1) * prod
    return 2 * (q**k + (-1) ** (k + 1)) * prod


def test_GLnq_order():
    assert GL(2, 3).order == get_GLnq_order(2, 3)


def test_SLnq_order():
    assert SL(2, 3).order == get_SLnq_order(2, 3)


def test_Onq_order():
    characteristic = 7
    degree = 1
    dimension = 2
    assert O(dimension, characteristic**degree).order == get_Onq_order(
        dimension, characteristic**degree
    )
    characteristic = 3
    degree = 2
    dimension = 2
    assert O(dimension, characteristic**degree).order == get_Onq_order(
        dimension, characteristic**degree
    )


def test_heisenberg_order():
    assert heisenberg_group(3, 5).order == 5**3
