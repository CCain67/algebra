"""Unit tests for group elements"""

import galois

from noetherpy.group_theory.group_elements import (
    CyclicGroupElement,
    Matrix,
    Permutation,
)

GF = galois.GF(2)


def test_cyclic_group_identity_elements() -> None:
    assert (
        CyclicGroupElement(7, 0).is_identity()
        == CyclicGroupElement(8, 0).is_identity()
        == CyclicGroupElement(9, 0).is_identity()
    )


def test_same_cyclic_group_elements_with_different_symbols_should_be_equal() -> None:
    a = CyclicGroupElement(8, 1, "a")
    x = CyclicGroupElement(8, 1, "x")
    assert a == x


def test_cyclic_group_element_to_permutation() -> None:
    assert CyclicGroupElement(7, 1).to_permutation() == Permutation(
        [1, 2, 3, 4, 5, 6, 0]
    )
    assert CyclicGroupElement(8, 2).to_permutation() == Permutation(
        [2, 3, 4, 5, 6, 7, 0, 1]
    )


def test_cyclic_group_element_to_matrix() -> None:
    assert CyclicGroupElement(4, 1).to_matrix() == Matrix(
        GF([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]]), 2, 1
    )
    assert CyclicGroupElement(6, 3).to_matrix() == Matrix(
        GF(
            [
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
            ]
        ),
        2,
        1,
    )


p = Permutation([1, 2, 0, 4, 3])
q = Permutation([1, 2, 3, 4, 0])
r = Permutation(list(range(5)))
s = Permutation([1, 0, 2, 4, 3])
t = Permutation([0, 3, 2, 1, 4])


def test_permutation_cycle_decomposition() -> None:
    assert p.cycle_decomposition == [[0, 1, 2], [3, 4]]
    assert q.cycle_decomposition == [[0, 1, 2, 3, 4]]
    assert r.cycle_decomposition == []
    assert s.cycle_decomposition == [[0, 1], [3, 4]]
    assert t.cycle_decomposition == [[1, 3]]


def test_permutation_cycle_notation() -> None:
    assert p.cycle_notation == "(3 4)(0 1 2)"
    assert q.cycle_notation == "(0 1 2 3 4)"
    assert r.cycle_notation == "1"
    assert s.cycle_notation == "(3 4)(0 1)"
    assert t.cycle_notation == "(1 3)"


def test_permutation_cycle_type() -> None:
    assert p.cycle_type == [2, 3]
    assert q.cycle_type == [5]
    assert r.cycle_type == []
    assert s.cycle_type == [2, 2]
    assert t.cycle_type == [2]


def test_permutation_order() -> None:
    assert p.order == 6
    assert q.order == 5
    assert r.order == 1
    assert s.order == 2
    assert t.order == 2


def test_permutation_sign() -> None:
    assert p.sign == -1
    assert q.sign == 1
    assert r.sign == 1
    assert s.sign == 1
    assert t.sign == -1


def test_permutation_is_cycle() -> None:
    assert p.is_identity() is False
    assert q.is_identity() is False
    assert r.is_identity() is True
    assert s.is_identity() is False
    assert t.is_identity() is False


def test_permutation_is_identity() -> None:
    assert p.is_cycle() is False
    assert q.is_cycle() is True
    assert r.is_cycle() is False
    assert s.is_cycle() is False
    assert t.is_cycle() is True


def test_permutation_is_transposition() -> None:
    assert p.is_transposition() is False
    assert q.is_transposition() is False
    assert r.is_transposition() is False
    assert s.is_transposition() is False
    assert t.is_transposition() is True


def test_permutation_to_matrix() -> None:
    assert p.to_matrix() == Matrix(
        GF(
            [
                [0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1],
                [0, 0, 0, 1, 0],
            ]
        ),
        2,
        1,
    )
    assert q.to_matrix() == Matrix(
        GF(
            [
                [0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1],
                [1, 0, 0, 0, 0],
            ]
        ),
        2,
        1,
    )
    assert r.to_matrix() == Matrix(
        GF(
            [
                [1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1],
            ]
        ),
        2,
        1,
    )
    assert s.to_matrix() == Matrix(
        GF(
            [
                [0, 1, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1],
                [0, 0, 0, 1, 0],
            ]
        ),
        2,
        1,
    )
    assert t.to_matrix() == Matrix(
        GF(
            [
                [1, 0, 0, 0, 0],
                [0, 0, 0, 1, 0],
                [0, 0, 1, 0, 0],
                [0, 1, 0, 0, 0],
                [0, 0, 0, 0, 1],
            ]
        ),
        2,
        1,
    )
