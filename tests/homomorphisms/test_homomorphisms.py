import pytest

from group_theory.base.group_elements import GroupElement
from group_theory.base.homomorphisms import Homomorphism
from group_theory.common_groups.permutation_groups import SymmetricGroup

G = SymmetricGroup.as_permutation_group(4)
g = G[2]

def conjugation(x: GroupElement) -> GroupElement:
    return g*x*(~g)
F = Homomorphism(G,conjugation,G)

def left_multiplication(x: GroupElement) -> GroupElement:
    return g*x
L = Homomorphism(G,left_multiplication,G)

def test_conjugation():
    assert F.validate_homomorphism()

def test_kernel():
    assert F.kernel.is_trivial()

def test_left_multiplication():
    assert ~L.validate_homomorphism()