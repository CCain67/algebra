"""
This module defines the Ring, Subring, and Ideal classes. 
"""

from __future__ import annotations

from ..group_theory.groups import Group
from .ring_elements import (
    GroupElementAdapter,
    RingElement,
)


class Ring:
    """Base class representing a finite ring.

    Args:
        elements (list[RingElement]): a list of elements for the ring.
    """

    # pylint: disable=too-many-instance-attributes
    # pylint: disable=too-many-public-methods

    def __init__(self, elements: list[RingElement]):
        self.elements = elements
        self.order = len(self.elements)
        self.additive_identity = self.get_additive_identity()
        self.multiplicative_identity = self.get_multiplicative_identity()

    def __repr__(self) -> str:
        repr_string = ""
        for ring_element in self:
            repr_string += repr(ring_element) + "\n"
        return repr_string

    def __eq__(self, other) -> bool:
        return set(self.elements) == set(other.elements)

    def __ne__(self, other) -> bool:
        return set(self.elements) != set(other.elements)

    def __lt__(self, other) -> bool:
        return set(self.elements).issubset(set(other.elements))

    def __le__(self, other) -> bool:
        return set(self.elements).issubset(set(other.elements))

    def __gt__(self, other) -> bool:
        return set(self.elements).issuperset(set(other.elements))

    def __ge__(self, other) -> bool:
        return set(self.elements).issuperset(set(other.elements))

    def __hash__(self) -> int:
        return hash(frozenset(self.elements))

    def __getitem__(self, key) -> RingElement:
        return self.elements[key]

    def __iter__(self):
        return iter(self.elements)

    def __len__(self) -> int:
        return len(self.elements)

    def get_additive_identity(self) -> RingElement:
        """Fetches the additive identity element of the ring

        Returns:
            RingElement: the additive identity element of the ring.
        """
        for element in self.elements:
            if element.is_additive_identity():
                self.identity = element
                return element
        return None

    def get_multiplicative_identity(self) -> RingElement:
        """Fetches the multiplicative identity element of the ring

        Returns:
            RingElement: the multiplicative identity element of the ring.
        """
        for element in self.elements:
            if element.is_multiplicative_identity():
                self.identity = element
                return element
        return None

    def group_of_units(self) -> Group:
        """Computes the group of units of the ring.

        Returns:
            Group: the group of invertible elements of the ring.
        """
        return Group([GroupElementAdapter(r) for r in self if r.is_invertible()])
