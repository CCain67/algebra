from group_theory.group_elements import (
    Permutation,
)
from group_theory.groups import Group
from group_theory.common_groups.residue_class_groups import cyclic_group

def klein_four_group(repr: str = 'residue') -> Group:
    if repr not in ['residue', 'permutation', 'matrix']:
        raise ValueError('repr must be one of: "residue", "permutation" or "matrix"')

    if repr=='residue':
        C_2 = cyclic_group(2)
        return C_2*C_2
    elif repr=='permutation' or repr=='matrix':
        return dihedral_group(2, repr)

def dihedral_group(N: int, repr: str = 'permutation') -> Group:
    if repr not in ['permutation', 'matrix']:
        raise ValueError('repr must be one of: "permutation" or "matrix"')
    
    cycle = {**{i:i+1 for i in range(1,N)}, N:1}
    flip = {1:1, **{2+i:N-i for i in range(N-1)}}

    if repr=='permutation':
        r = Permutation(cycle)
        s = Permutation(flip)
    elif repr=='matrix':
        r = Permutation(cycle).to_matrix()
        s = Permutation(flip).to_matrix()

    D = [(r**i)*(s**j) for i in range(N) for j in range(2)]
    D_n = Group(D)
    D_n.canonical_generators = [r,s]
    return D_n

def quasidihedral_group(N: int, repr: str) -> Group:
    return NotImplemented

def dicyclic_group(N: int, repr: str) -> Group:
    return NotImplemented