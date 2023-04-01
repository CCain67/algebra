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

class CyclicGroup(Group):
    '''
    Basic abstract factory for cyclic groups. 
    '''
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)
    
    @classmethod
    def as_residue_class_group(cls, modulus: int):
        g = AdditiveResidueClass(1, modulus)
        C = cls([g**j for j in range(modulus)])
        C.canonical_generators = [g]
        return C
    
    @classmethod
    def as_permutation_group(cls, modulus: int):
        generator = {i:i+1 for i in range(1,modulus)}
        generator[modulus]=1
        g = Permutation(generator)
        
        C = cls([g**j for j in range(modulus)])
        C.canonical_generators = [g]
        return C

    @classmethod
    def as_matrix_group(cls, modulus: int):
        generator = {i:i+1 for i in range(1,modulus)}
        generator[modulus]=1
        M = numpy.zeros((modulus, modulus))
        for k in generator.keys():
            M[k-1,generator[k]-1]=1
        g = Matrix(M)

        C = cls([g**j for j in range(modulus)])
        C.canonical_generators = [g]
        return C
    
    
class GroupOfUnits(Group):
    '''
    Basic abstract factory for cyclic groups. 
    '''
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)
    
    @classmethod
    def as_residue_class_group(cls, modulus: int):
        return cls([MultiplicativeResidueClass(i, modulus) for i in range(modulus) if math.gcd(i,modulus)==1])