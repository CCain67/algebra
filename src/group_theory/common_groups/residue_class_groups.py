import math
import numpy

from group_theory.base.group_elements import (
    AdditiveResidueClass,
    GroupElement,
    Matrix,
    MultiplicativeResidueClass,
    Permutation,
    ResidueClass,
)
from group_theory.base.groups import Group

def cyclic_group(N: int, repr: str = 'residue') -> Group:
    if repr not in ['residue', 'permutation', 'matrix']:
        raise ValueError('repr must be one of: "residue", "permutation" or "matrix"')
    
    if repr=='residue':
        g = AdditiveResidueClass(1, N)
        C = Group([g**j for j in range(N)])
        C.canonical_generators = [g]
        return C
    
    elif repr=='permutation':
        generator = {i:i+1 for i in range(1,N)}
        generator[N]=1
        g = Permutation(generator)
        
        C = Group([g**j for j in range(N)])
        C.canonical_generators = [g]
        return C
    
    elif repr=='matrix':
        generator = {i:i+1 for i in range(1,N)}
        generator[N]=1
        M = numpy.zeros((N, N))
        for k in generator.keys():
            M[k-1,generator[k]-1]=1
        g = Matrix(M)

        C = Group([g**j for j in range(N)])
        C.canonical_generators = [g]
        return C