import itertools

from group_theory.base.group_elements import (
    GroupElement,
    Permutation,
)
from group_theory.base.groups import Group

'''
below are basic abstract factories for the symmetric and alternating groups. 
'''

class SymmetricGroup(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)
    
    @classmethod
    def as_permutation_group(cls, N: int):
        if N<1:
            raise ValueError('the symmetric group on N letters requires N to be at least 1')
        
        P = list(itertools.permutations(range(1,N+1)))
        S_N = cls([Permutation({i+1:p[i] for i in range(N)}) for p in P])
        
        t = {1:2, 2:1}
        if N==2:
            S_N.canonical_generators = [Permutation(t)]
        elif N>2:
            for i in range(3,N+1):
                t[i]=i
            c = {i:i+1 for i in range(1,N)}
            c[N]=1
            S_N.canonical_generators = [Permutation(c), Permutation(t)]
        return S_N
    
    @classmethod
    def as_matrix_group(cls, N: int):
        if N<1:
            raise ValueError('the symmetric group on N letters requires N to be at least 1')
        P = list(itertools.permutations(range(1,N+1)))
        S_N = cls([Permutation({i+1:p[i] for i in range(N)}).to_matrix() for p in P])
        
        t = {1:2, 2:1}
        if N==2:
            S_N.canonical_generators = [Permutation(t).to_matrix()]
        elif N>2:
            for i in range(3,N+1):
                t[i]=i
            c = {i:i+1 for i in range(1,N)}
            c[N]=1
            S_N.canonical_generators = [Permutation(c).to_matrix(), Permutation(t).to_matrix()]
        return S_N
    
class AlternatingGroup(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)
    
    @classmethod
    def as_permutation_group(cls, N: int):
        P = list(itertools.permutations(range(1,N+1)))
        Sym = [Permutation({i+1:p[i] for i in range(N)}) for p in P]
        return cls([x for x in Sym if x.sign==1])