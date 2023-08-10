"""Unit tests for Group, Subgroup, and Coset classes"""

from noetherpy.group_theory.group_elements import GroupElement, Permutation
from noetherpy.group_theory.groups import Coset, Group, Subgroup

from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group


def test_trivial_group_properties() -> None:
    class TrivialGroupElement(GroupElement):
        """Abstract trivial element"""

        def __mul__(self, other):
            return self

        def __invert__(self):
            return self

        def __eq__(self, other):
            if isinstance(other, TrivialGroupElement):
                return True
            return False

        def __hash__(self):
            return hash("trivial group element")

        def get_order(self):
            return 1

        def is_identity(self):
            return True

    trivial_group = Group([TrivialGroupElement()])
    trivial_subgroup = Subgroup([TrivialGroupElement()], trivial_group)

    # group properties
    assert trivial_group.is_trivial()
    assert trivial_group.is_abelian()
    assert trivial_group.is_perfect()

    # subgroup properties
    assert trivial_subgroup.validate_inclusion()
    assert trivial_subgroup.validate_subgroup()
    assert trivial_subgroup.is_normal
    assert (
        trivial_group.subgroup_generated_by([TrivialGroupElement()]) == trivial_subgroup
    )

    # properties of operations
    assert (trivial_group * trivial_group).is_trivial()
    assert (trivial_group / trivial_subgroup).is_trivial()


def test_cosets_should_partition_the_group() -> None:
    G = symmetric_group(3)
    H = G.subgroup_generated_by([Permutation([1, 0, 2])])
    cosets = G % H

    assert not H.is_normal
    assert len(cosets) == G.order / H.order
    assert set.intersection(*[set(X.elements) for X in cosets]) == set()
    assert set.union(*[set(X.elements) for X in cosets]) == set(G.elements)
