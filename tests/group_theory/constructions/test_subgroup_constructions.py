"""Unit tests for subgroup constructions"""
from noetherpy.group_theory.groups import Subgroup

from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group

from noetherpy.group_theory.constructions.subgroup_constructions import (
    centralizer,
    conjugate_subgroup,
    normal_closure,
    normal_core,
    normalizer,
)


def test_normalizer():
    G = symmetric_group(5)
    S = G.subgroup_generated_by([G.canonical_generators[0]])
    N = normalizer(S, G)
    C = centralizer(S, G)
    S = Subgroup(S.elements, N)
    assert S.is_normal
    assert C < N
