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

def linear_combinations(field, vector_list: list):
    number_of_vectors = len(vector_list)
    all_scalar_combinations = list(product(*[[field(i) for i in range(field.order)]]*number_of_vectors))
    summand_list = []
    for combination in all_scalar_combinations:
        summand_list += [[c*x for c,x in list(zip(combination,vector_list))]]
    linear_combinations = [reduce(lambda x,y: x+y, summand) for summand in summand_list]
    return linear_combinations

def remove_linear_combos(field, starting_vectors, vector_list):
    linear_combos = linear_combinations(field,vector_list)
    return [z for z in starting_vectors if all(any(z!=y) for y in linear_combos)]

class GeneralLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        '''
        this is essentially an implementation of the standard proof of the order of GL(n,q): an invertible matrix is produced row by row,
        removing all linear combinations of the previous rows before selecting the next one, until all rows are filled. 
        '''
        GF = galois.GF(characteristic,degree, repr="poly")
        starting_vector_arrays = list(product(*[range(characteristic**degree)]*dimension))[1:]
        starting_vectors = [GF(a) for a in starting_vector_arrays]
        GL = []
        for v in starting_vectors:
            k = 0
            v_lists = [[v]]
            while k<dimension-1:
                v_lists = [l+[w] for l in v_lists for w in remove_linear_combos(GF,starting_vectors,l)]
                k+=1
            GL += [Matrix(GF(m), characteristic, degree) for m in v_lists]
        GLnq = cls(GL)
        GLnq.ground_field = GF
        return GLnq
    

class ProjectiveGeneralLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, dimension: int, characteristic: int, degree: int):
        '''
        all projective matrix groups are quotients by the center: the center happens to always be scalar multiples of the identity, 
        which is faster to define than calling MatrixGroup.center
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
        starting_vector_arrays = list(product(*[range(characteristic**degree)]*dimension))[1:]
        starting_vectors = [GF(a) for a in starting_vector_arrays]
        GL = []
        for v in starting_vectors:
            k = 0
            v_lists = [[v]]
            while k<dimension-1:
                v_lists = [l+[w] for l in v_lists for w in remove_linear_combos(GF,starting_vectors,l)]
                k+=1
            GL += [GF(m) for m in v_lists]
        SL = [Matrix(M, characteristic, degree) for M in GL if np.linalg.det(M)==1]
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
        GF = galois.GF(characteristic,degree, repr="poly")
        starting_vector_arrays = list(product(*[range(characteristic**degree)]*dimension))[1:]
        starting_vectors = [GF(a) for a in starting_vector_arrays]
        GL = []
        for v in starting_vectors:
            k = 0
            v_lists = [[v]]
            while k<dimension-1:
                v_lists = [l+[w] for l in v_lists for w in remove_linear_combos(GF,starting_vectors,l)]
                k+=1
            GL += [GF(m) for m in v_lists]
        O = [Matrix(M, characteristic, degree) for M in GL if (M*M.transpose()==GF.Identity(dimension)).all()]
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
        starting_vector_arrays = list(product(*[range(characteristic**degree)]*dimension))[1:]
        starting_vectors = [GF(a) for a in starting_vector_arrays]
        GL = []
        for v in starting_vectors:
            k = 0
            v_lists = [[v]]
            while k<dimension-1:
                v_lists = [l+[w] for l in v_lists for w in remove_linear_combos(GF,starting_vectors,l)]
                k+=1
            GL += [GF(m) for m in v_lists]
        O = [M for M in GL if (M*M.transpose()==GF.Identity(dimension)).all()]
        SO = [Matrix(M, characteristic, degree) for M in O if np.linalg.det(M)==1]
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