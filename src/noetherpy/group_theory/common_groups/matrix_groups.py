"""This module defines functions for constructing various matrix groups over finite fields."""

from typing import Type
from itertools import product
from functools import reduce

import numpy as np
import galois

from ..group_elements import (
    Matrix,
)
from ..groups import (
    Group,
    LinearGroup,
    Subgroup,
)


def _linear_combinations(GF: Type[galois.FieldArray], vector_list: list) -> list:
    """Constructs all linear combinations of the vectors given over the finite field.

    Args:
        GF (Type[galois.FieldArray]): ground field
        vector_list (list): list of vectors to form linear combinations from

    Returns:
        list: list of all possible linear combinations of the vectors provided
    """
    number_of_vectors = len(vector_list)
    all_scalar_combinations = product(
        *[[GF(i) for i in range(GF.order)]] * number_of_vectors
    )
    summand_list = []
    for combination in all_scalar_combinations:
        summand_list += [[c * x for c, x in zip(combination, vector_list)]]
    return [reduce(lambda x, y: x + y, summand) for summand in summand_list]


def _remove_linear_combos(
    GF: Type[galois.FieldArray], starting_vectors: list, vector_list: list
) -> list:
    """Filters out linear combinations of vectors from list of starting vectors

    Args:
        GF (Type[galois.FieldArray]): ground field
        starting_vectors (list): list of vectors to filter
        vector_list (list): list of vectors whose linear combinations need to be removed

    Returns:
        list: starting_vector list with linear combinations of vectors in vector_list removed
    """
    linear_combos = _linear_combinations(GF, vector_list)
    return [z for z in starting_vectors if all(any(z != y) for y in linear_combos)]


def _get_general_linear_matrices(
    GF: Type[galois.FieldArray], dimension: int
) -> list[Matrix]:
    """
    this is essentially an implementation of the standard proof of the order of GL(n,q):
    an invertible matrix is produced row by row, removing all linear combinations of the
    previous rows before selecting the next one, until all rows are filled.
    """
    characteristic = GF.characteristic
    degree = GF.degree
    starting_vector_entries = list(
        product(*[range(characteristic**degree)] * dimension)
    )[
        1:
    ]  # the 0th tuple is always (0,0,...,0)
    starting_vectors = [GF(a) for a in starting_vector_entries]
    general_linear_matrices = []
    for vector in starting_vectors:
        k = 0
        v_lists = [[vector]]
        while k < dimension - 1:
            v_lists = [
                l + [w]
                for l in v_lists
                for w in _remove_linear_combos(GF, starting_vectors, l)
            ]
            k += 1
        general_linear_matrices += [
            Matrix(GF(m), characteristic, degree) for m in v_lists
        ]
    return general_linear_matrices


def _filter_determinant_one(matrix_list: list[Matrix]) -> list[Matrix]:
    return [M for M in matrix_list if np.linalg.det(M.matrix) == 1]


def _filter_orthogonal(
    matrix_list: list[Matrix], GF: Type[galois.FieldArray], dimension: int
) -> list[Matrix]:
    return [
        M
        for M in matrix_list
        if (M.matrix @ M.matrix.transpose() == GF.Identity(dimension)).all()
    ]


def _get_central_elements(
    matrix_list: list[Matrix], GF: Type[galois.FieldArray], dimension: int
) -> list[Matrix]:
    return [
        Matrix(c * GF.Identity(dimension), GF.characteristic, GF.degree)
        for c in GF.elements
        if Matrix(c * GF.Identity(dimension), GF.characteristic, GF.degree)
        in matrix_list
    ]


def GL(n: int, q: int) -> LinearGroup:
    """Constructs the general linear group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the general linear of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    return LinearGroup(_get_general_linear_matrices(GF, n), GF)


def SL(n: int, q: int) -> LinearGroup:
    """Constructs the special linear group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the special linear of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    general_linear_matrices = _get_general_linear_matrices(GF, n)
    return LinearGroup(_filter_determinant_one(general_linear_matrices), GF)


def PGL(n: int, q: int) -> LinearGroup:
    """Constructs the projective general linear group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the projective general linear group of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    general_linear_group = Group(_get_general_linear_matrices(GF, n))
    gl_center = Subgroup(
        _get_central_elements(general_linear_group.elements, GF, n),
        general_linear_group,
    )
    return LinearGroup((general_linear_group / gl_center).elements, GF)


def PSL(n: int, q: int) -> LinearGroup:
    """Constructs the projective special linear group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the projective special linear group of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    general_linear_matrices = _get_general_linear_matrices(GF, n)
    special_linear_group = Group(_filter_determinant_one(general_linear_matrices))
    sl_center = Subgroup(
        _get_central_elements(special_linear_group.elements, GF, n),
        special_linear_group,
    )
    return LinearGroup((special_linear_group / sl_center).elements, GF)


def O(n: int, q: int) -> LinearGroup:
    """Constructs the orthogonal group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the orthogonal group of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    general_linear_matrices = _get_general_linear_matrices(GF, n)
    return LinearGroup(_filter_orthogonal(general_linear_matrices, GF, n), GF)


def SO(n: int, q: int) -> LinearGroup:
    """Constructs the special orthogonal group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the special orthogonal group of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    general_linear = _get_general_linear_matrices(GF, n)
    special_linear = _filter_determinant_one(general_linear)
    return LinearGroup(_filter_orthogonal(special_linear, GF, n), GF)


def PO(n: int, q: int) -> LinearGroup:
    """Constructs the projective orthogonal group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the projective orthogonal group of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    general_linear = _get_general_linear_matrices(GF, n)
    orthogonal_group = Group(_filter_orthogonal(general_linear, GF, n))
    orthogonal_center = Subgroup(
        _get_central_elements(orthogonal_group.elements, GF, n), orthogonal_group
    )
    return LinearGroup((orthogonal_group / orthogonal_center).elements, GF)


def PSO(n: int, q: int) -> LinearGroup:
    """Constructs the projective special orthogonal group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the projective special orthogonal group of dimension n over
        the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    general_linear_matrices = _get_general_linear_matrices(GF, n)
    special_linear_matrices = _filter_determinant_one(general_linear_matrices)
    special_orthogonal_group = Group(_filter_orthogonal(special_linear_matrices, GF, n))
    so_center = Subgroup(
        _get_central_elements(special_orthogonal_group.elements, GF, n),
        special_orthogonal_group,
    )
    return LinearGroup((SO / so_center).elements, GF)


def _tri(n: int) -> int:
    """Returns the nth triangular number"""
    return sum(range(1, n + 1))


def _tuple_to_matrix(
    dimension: int, GF: Type[galois.FieldArray], element_tuple: tuple
) -> galois.FieldArray:
    """
    Creates an element of th Heisenberg group from a vector of finite field elements.

    Note: here, we need len(t)==tri(dimension)
    """
    heisenberg_matrix = GF.Identity(dimension)
    counter = 0
    for pair in [
        pairs
        for pairs in product(range(dimension), range(dimension))
        if pairs[1] > pairs[0]
    ]:
        heisenberg_matrix[pair[0], pair[1]] = element_tuple[counter]
        counter += 1
    return heisenberg_matrix


def heisenberg_group(n: int, q: int) -> LinearGroup:
    """Constructs the Heisenberg group of dimension n over a finite field.

    Args:
        n (int): the dimension of the matrices in the group
        q (int): the prime power order of the ground field

    Returns:
        LinearGroup: the Heisenberg group of dimension n over the finite field F_q
    """
    GF = galois.GF(q, repr="poly")
    all_upper_entries = [GF(e) for e in list(product(*[range(GF.order)] * _tri(n - 1)))]
    heisenberg_matrices = [_tuple_to_matrix(n, GF, t) for t in all_upper_entries]
    return LinearGroup(
        [Matrix(M, GF.characteristic, GF.degree) for M in heisenberg_matrices], GF
    )
