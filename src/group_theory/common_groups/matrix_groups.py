import numpy as np
import sympy
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
        return cls(GL)
    

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
        starting_vectors = set(product(*[range(p)]*n))-{tuple([0]*n)}
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
        starting_vectors = set(product(*[range(p)]*n))-{tuple([0]*n)}
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
    

def tri(n: int) -> int:
    return sum(range(1,n+1))

def tuple_to_matrix(n: int, t: tuple) -> sympy.Matrix:
    '''
    here, we need len(t)==tri(n)
    '''
    m = sympy.Matrix.eye(n)
    counter = 0
    for pair in [pairs for pairs in product(range(n),range(n)) if pairs[1]>pairs[0]]:
        m[pair[0],pair[1]] =  t[counter]
        counter+=1
    return m

class HeisenbergGroup(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def over_finite_field(cls, n: int, p: int):
        all_upper_entries = list(product(*[range(p)]*tri(n-1)))
        H = [tuple_to_matrix(n,t) for t in all_upper_entries]
        return cls([Matrix(M,p) for M in H])