import numpy as np
import sympy

import itertools
from functools import reduce

from group_theory.base.group_elements import (
    GroupElement,
    Matrix,
)
from group_theory.base.groups import (
    Group,
    Subgroup,
)

def scalar_multiplication(p: int, c: int, vec: tuple):
    return tuple((c*x)%p for x in vec)

def vector_addition(p: int, vec_1: tuple, vec_2: tuple):
    return tuple((x+y)%p for x,y in zip(vec_1, vec_2))

def linear_combinations(p: int, vector_list: list):
    number_of_vectors = len(vector_list)
    all_scalar_combinations = list(itertools.product(*[range(p)]*number_of_vectors))
    summand_list = []
    for combination in all_scalar_combinations:
        summand_list += [[scalar_multiplication(p,*x) for x in list(zip(combination,vector_list))]]
    linear_combinations = {reduce(lambda x,y: vector_addition(p,x,y), summand) for summand in summand_list}
    return linear_combinations

class GeneralLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, n: int, p: int):
        starting_vectors = set(itertools.product(*[range(p)]*n))-{tuple([0]*n)}
        GL = []
        for v in starting_vectors:
            k = 0
            v_lists = [[v]]
            while k<n-1:
                v_lists = [l+[w] for l in v_lists for w in starting_vectors-linear_combinations(p,l)]
                k+=1
            GL += [sympy.Matrix(m) for m in v_lists]
        return cls([Matrix(M,p) for M in GL])
    

class ProjectiveGeneralLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, n: int, p: int):
        GL = GeneralLinear.over_finite_field(n=n,p=p)
        z = [M for M in GL if M.matrix.is_diagonal() and tuple(M.matrix.diagonal()) in {tuple([i]*n) for i in range(1,p)}]
        Z = Subgroup(z,GL)
        return GL/Z
    
class SpecialLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, n: int, p: int):
        starting_vectors = set(itertools.product(*[range(p)]*n))-{tuple([0]*n)}
        GL = []
        for v in starting_vectors:
            k = 0
            v_lists = [[v]]
            while k<n-1:
                v_lists = [l+[w] for l in v_lists for w in starting_vectors-linear_combinations(p,l)]
                k+=1
            GL += [sympy.Matrix(m) for m in v_lists]
        return cls([Matrix(M,p) for M in GL if M.det()%p==1])
    

class ProjectiveSpecialLinear(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, n: int, p: int):
        SL = SpecialLinear.over_finite_field(n=n,p=p)
        z = [M for M in SL if M.matrix.is_diagonal() and tuple(M.matrix.diagonal()) in {tuple([i]*n) for i in range(1,p)}]
        Z = Subgroup(z,SL)
        return SL/Z
    
class Orthogonal(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, n: int, p: int):
        starting_vectors = set(itertools.product(*[range(p)]*n))-{tuple([0]*n)}
        GL = []
        for v in starting_vectors:
            k = 0
            v_lists = [[v]]
            while k<n-1:
                v_lists = [l+[w] for l in v_lists for w in starting_vectors-linear_combinations(p,l)]
                k+=1
            GL += [sympy.Matrix(m) for m in v_lists]
        return cls([Matrix(M,p) for M in GL if M*M.transpose()%p==sympy.eye(n)])
    

class ProjectiveOrthogonal(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, n: int, p: int):
        O = Orthogonal.over_finite_field(n=n,p=p)
        z = [M for M in O if M.matrix.is_diagonal() and tuple(M.matrix.diagonal()) in {tuple([i]*n) for i in range(1,p)}]
        Z = Subgroup(z,O)
        return O/Z