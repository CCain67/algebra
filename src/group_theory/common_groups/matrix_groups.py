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

def _get_GLnq_matrices(GF: Type[galois.FieldArray], dimension: int) -> list:
    '''
    this is essentially an implementation of the standard proof of the order of GL(n,q): an invertible matrix is produced row by row,
    removing all linear combinations of the previous rows before selecting the next one, until all rows are filled. 
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

class GeneralLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        GF = galois.GF(characteristic,degree, repr="poly")
        GLnq = cls(_get_GLnq_matrices(GF, dimension))
        GLnq.ground_field = GF
        return GLnq
    

class ProjectiveGeneralLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        '''
        all projective matrix groups are quotients by the center: the center happens to always be scalar multiples of the identity, 
        which is faster to define than calling Group.center
        '''
        GL = GeneralLinear.over_finite_field(dimension, characteristic, degree)
        Z = Subgroup([Matrix(c*GL.ground_field.Identity(dimension),characteristic,degree) for c in GL.ground_field.elements],GL)
        return GL/Z
    
class SpecialLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        GF = galois.GF(characteristic,degree, repr="poly")
        GL = _get_GLnq_matrices(GF, dimension)
        SL = [M for M in GL if np.linalg.det(M.matrix)==1]
        SLnq = cls(SL)
        SLnq.ground_field = GF
        return SLnq
    

class ProjectiveSpecialLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        SL = SpecialLinear.over_finite_field(dimension, characteristic, degree)
        Z = Subgroup([Matrix(c*SL.ground_field.Identity(dimension),characteristic,degree) for c in SL.ground_field.elements],SL)
        return SL/Z
    
class Orthogonal(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        
        GF = galois.GF(characteristic, degree, repr="poly")
        GL = _get_GLnq_matrices(GF, dimension)
        O = [M for M in GL if (M.matrix@M.matrix.transpose()==GF.Identity(dimension)).all()]
        Onq = cls(O)
        Onq.ground_field = GF
        return Onq
    
class ProjectiveOrthogonal(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        O = Orthogonal.over_finite_field(dimension, characteristic, degree)
        Z = Subgroup([Matrix(c*O.ground_field.Identity(dimension),characteristic,degree) for c in O.ground_field.elements],O)
        return O/Z
    
class SpecialOrthogonal(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        GF = galois.GF(characteristic,degree, repr="poly")
        GL = _get_GLnq_matrices(GF, dimension)
        O = [M for M in GL if (M.matrix@M.matrix.transpose()==GF.Identity(dimension)).all()]
        SO = [M for M in O if np.linalg.det(M.matrix)==1]
        SOnq = cls(SO)
        SOnq.ground_field = GF
        return SOnq
    
class ProjectiveSpecialOrthogonal(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        SO = SpecialOrthogonal.over_finite_field(dimension, characteristic, degree)
        Z = Subgroup([Matrix(c*SO.ground_field.Identity(dimension),characteristic,degree) for c in SO.ground_field.elements],SO)
        return SO/Z

def tri(n: int) -> int:
    return sum(range(1,n+1))

def tuple_to_matrix(dimension: int, GF: Type[galois.FieldArray], t: tuple) -> galois.FieldArray:
    '''
    here, we need len(t)==tri(dimension)
    '''
    M = GF.Identity(dimension)
    counter = 0
    for pair in [pairs for pairs in product(range(dimension),range(dimension)) if pairs[1]>pairs[0]]:
        M[pair[0],pair[1]] =  t[counter]
        counter+=1
    return M

class HeisenbergGroup(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        GF = galois.GF(characteristic,degree, repr="poly")
        all_upper_entries = [GF(e) for e in list(product(*[range(GF.order)]*tri(dimension-1)))]
        H = [tuple_to_matrix(dimension, GF, t) for t in all_upper_entries]
        return cls([Matrix(M, characteristic, degree) for M in H])