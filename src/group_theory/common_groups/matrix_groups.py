from typing import Type

import numpy as np
import galois

from itertools import product
from functools import reduce

from group_theory.base.group_elements import (
    GroupElement,
    Matrix,
)
from group_theory.base.groups import (
    Group,
    Subgroup,
)

def _linear_combinations(GF: Type[galois.FieldArray], vector_list: list) -> list:
    number_of_vectors = len(vector_list)
    all_scalar_combinations = product(*[[GF(i) for i in range(GF.order)]]*number_of_vectors)
    summand_list = []
    for combination in all_scalar_combinations:
        summand_list += [[c*x for c,x in zip(combination,vector_list)]]
    return [reduce(lambda x,y: x+y, summand) for summand in summand_list]

def _remove_linear_combos(GF: Type[galois.FieldArray], starting_vectors: list, vector_list: list) -> list:
    linear_combos = _linear_combinations(GF,vector_list)
    return [z for z in starting_vectors if all(any(z!=y) for y in linear_combos)]

def _get_GLnq_matrices(GF: Type[galois.FieldArray], dimension: int) -> list[Matrix]:
    '''
    this is essentially an implementation of the standard proof of the order of GL(n,q): 
    an invertible matrix is produced row by row, removing all linear combinations of the 
    previous rows before selecting the next one, until all rows are filled. 
    '''
    characteristic = GF.characteristic
    degree = GF.degree
    starting_vector_entries = list(product(*[range(characteristic**degree)]*dimension))[1:] # the 0th tuple is always (0,0,...,0)
    starting_vectors = [GF(a) for a in starting_vector_entries]
    GL = []
    for v in starting_vectors:
        k = 0
        v_lists = [[v]]
        while k<dimension-1:
            v_lists = [l+[w] for l in v_lists for w in _remove_linear_combos(GF,starting_vectors,l)]
            k+=1
        GL += [Matrix(GF(m), characteristic, degree) for m in v_lists]
    return GL

def _filter_determinant_one(matrix_list: list[Matrix]) -> list[Matrix]:
    return [M for M in matrix_list if np.linalg.det(M.matrix)==1]

def _filter_orthogonal(matrix_list: list[Matrix], GF: Type[galois.FieldArray], dimension: int) -> list[Matrix]:
    return [M for M in matrix_list if (M.matrix@M.matrix.transpose()==GF.Identity(dimension)).all()]

def _get_central_elements(matrix_list: list[Matrix], GF: Type[galois.FieldArray], dimension: int) -> list[Matrix]:
    return [
        Matrix(c*GF.Identity(dimension),GF.characteristic, GF.degree) for c in GF.elements 
        if Matrix(c*GF.Identity(dimension),GF.characteristic,GF.degree) in matrix_list
    ]


class LinearGroup(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)


def GL(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GLnq = LinearGroup(_get_GLnq_matrices(GF,n))
    GLnq.ground_field = GF
    return GLnq

def SL(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GL = _get_GLnq_matrices(GF,n)
    SLnq = LinearGroup(_filter_determinant_one(GL))
    SLnq.ground_field = GF
    return SLnq

def PGL(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GL = Group(_get_GLnq_matrices(GF,n))
    Z = Subgroup(_get_central_elements(GL.elements, GF, n), GL)
    PGLnq = LinearGroup(
        (GL/Z).elements
    )
    PGLnq.ground_field = GF
    return PGLnq

def PSL(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GL = _get_GLnq_matrices(GF,n)
    SL = Group(_filter_determinant_one(GL))
    Z = Subgroup(_get_central_elements(SL.elements, GF, n), SL)
    PSLnq = LinearGroup(
        (SL/Z).elements
    )
    PSLnq.ground_field = GF
    return PSLnq

def O(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GL = _get_GLnq_matrices(GF,n)
    Onq = LinearGroup(_filter_orthogonal(GL, GF, n))
    Onq.ground_field = GF
    return Onq

def SO(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GL = _get_GLnq_matrices(GF,n)
    SL = _filter_determinant_one(GL)
    SOnq = LinearGroup(_filter_orthogonal(SL, GF, n))
    SOnq.ground_field = GF
    return SOnq

def PO(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GL = _get_GLnq_matrices(GF,n)
    O = Group(_filter_orthogonal(GL, GF, n))
    Z = Subgroup(_get_central_elements(O.elements, GF, n), O)
    POnq = LinearGroup(
        (O/Z).elements
    )
    POnq.ground_field = GF
    return POnq

def PSO(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr='poly')
    GL = _get_GLnq_matrices(GF,n)
    SL = _filter_determinant_one(GL)
    SO = Group(_filter_orthogonal(SL, GF, n))
    Z = Subgroup(_get_central_elements(SO.elements, GF, n), SO)
    PSOnq = LinearGroup(
        (SO/Z).elements
    )
    PSOnq.ground_field = GF
    return PSOnq


def _tri(n: int) -> int:
    return sum(range(1,n+1))

def _tuple_to_matrix(dimension: int, GF: Type[galois.FieldArray], t: tuple) -> galois.FieldArray:
    '''
    here, we need len(t)==tri(dimension)
    '''
    M = GF.Identity(dimension)
    counter = 0
    for pair in [pairs for pairs in product(range(dimension),range(dimension)) if pairs[1]>pairs[0]]:
        M[pair[0],pair[1]] =  t[counter]
        counter+=1
    return M

def HeisenbergGroup(n: int, q: int) -> LinearGroup:
    GF = galois.GF(q, repr="poly")
    all_upper_entries = [GF(e) for e in list(product(*[range(GF.order)]*_tri(n-1)))]
    H = [_tuple_to_matrix(n, GF, t) for t in all_upper_entries]
    Hnq = LinearGroup([Matrix(M, GF.characteristic, GF.degree) for M in H])
    Hnq.ground_field = GF
    return Hnq