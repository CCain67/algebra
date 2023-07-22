"""
This module defines several types of products that can be formed between two groups:
- semidirect products
- central products
- wreath products
"""

from itertools import product
from typing import Callable

from group_theory.groups import (
    Group,
)
from group_theory.group_elements import GroupElement


class SemidirectProductElement(GroupElement):
    """
    Here, we are viewing the 'twisting' homomorphism twist:G -> Aut(H) as a map
    twist:G x H -> H which restricts to an automorphism of H for each g in G.
    """

    def __init__(
        self,
        elements: tuple[GroupElement],
        twist: Callable[[GroupElement, GroupElement], GroupElement],
    ) -> None:
        self.elements = elements
        self.twist = twist
        self.num_elements = 2

        # properties
        self._order = None

    def __repr__(self):
        return str(self.elements)

    def __hash__(self):
        return hash((self.elements, self.twist))

    def __eq__(self, other):
        return self.elements == other.elements

    def __ne__(self, other):
        return self.elements != other.elements

    def __mul__(self, other):
        if self.twist != other.twist:
            raise ValueError("the twisting homomorphisms do not match")
        element_product = (
            self.elements[0] * self.twist(self.elements[1], other.elements[0]),
            self.elements[1] * other.elements[1],
        )
        return SemidirectProductElement(element_product, self.twist)

    def __invert__(self):
        return SemidirectProductElement(
            (self.twist(~self.elements[1], ~self.elements[0]), ~self.elements[1]),
            self.twist,
        )

    def __getitem__(self, key):
        return self.elements[key]

    def __iter__(self):
        return iter(self.elements)

    def get_order(self) -> int:
        prod = self
        i = 1
        while not prod.is_identity():
            prod *= self
            i += 1
        return i

    @property
    def order(self) -> int:
        """Fetches the order of the element"""
        if self._order is None:
            self._order = self.get_order()
            return self._order
        return self._order

    def is_identity(self):
        return self.elements[0].is_identity() and self.elements[1].is_identity()


def semidirect_product(
    fiber: Group, twist: Callable[[GroupElement, GroupElement], GroupElement], base: Group
) -> Group:
    """Constructs the semidirect product of two groups given a group action.

    Args:
        - fiber (Group): the group being acted on
        - twist (Callable): the action of the base group on the fiber group
        - base (Group): the group which is performing the action.

    Terminology: for a group G with action f:G x H -> H, we will refer to G as the "base" group 
        and H as the "fiber" group, in an analogy to fiber bundles.

    Returns:
        Group: the semidirect product of the two groups provided via the action provided.
    """
    semidirect_product_elements = list(
        SemidirectProductElement(p, twist)
        for p in product(fiber.elements, base.elements)
    )
    return Group(semidirect_product_elements)
