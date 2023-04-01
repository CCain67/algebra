from group_theory.base.group_elements import (
    GroupElement,
    Matrix,
    Permutation,
    ResidueClass,
)
from group_theory.base.groups import Group
from group_theory.common_groups.residue_class_groups import CyclicGroup

class KleinFourGroup(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def as_residue_class_group(cls):
        C_2 = CyclicGroup.as_residue_class_group(2)
        return C_2*C_2
    
    @classmethod
    def as_permutation_group(cls):
        return DihedralGroup.as_permutation_group(2)

class QuaternionGroup(Group):
    pass

class DihedralGroup(Group):
    def __init__(self, elements: list[GroupElement]):
        super().__init__(elements)

    @classmethod
    def as_permutation_group(cls, n: int):
        cycle = {i:i+1 for i in range(1,n)}
        cycle[n]=1
        flip = {1:1}
        flip = {1:1}
        for i in range(n-1):
            flip[2+i]=n-i
        r = Permutation(cycle)
        s = Permutation(flip)

        D = [(r**i)*(s**j) for i in range(n) for j in range(2)]
        D_n = cls(D)
        D_n.canonical_generators = [r,s]
        return D_n


class QuasiDihedralGroup(Group):
    pass

class DicyclicGroup(Group):
    pass