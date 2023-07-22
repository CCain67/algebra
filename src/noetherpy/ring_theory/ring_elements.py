"""This module defines several various classes of ring elements."""

from __future__ import annotations
from abc import ABC, abstractmethod
from functools import reduce

import galois

from ..group_theory.group_elements import GroupElement


class RingElement(ABC):
    """Base class representing elements of arbitrary finite rings."""

    @abstractmethod
    def additive_order(self):
        """abstract method for fetching the additive order of the element."""

    @abstractmethod
    def multiplicative_order(self):
        """abstract method for fetching the additive order of the element."""

    @abstractmethod
    def is_additive_identity(self):
        """abstract method for determining whether or not the RingElement is
        the additive identity element of the ring.
        """

    @abstractmethod
    def is_multiplicative_identity(self):
        """abstract method for determining whether or not the RingElement is
        the multiplicative identity element of the ring.
        """

    @abstractmethod
    def is_invertible(self):
        """abstract method for determining whether or not the RingElement
        has a multiplicative inverse."""


class GroupElementAdapter(GroupElement):
    """An adapter which converts a RingElement to a GroupElement"""

    def __init__(self, ring_element: RingElement) -> None:
        self.ring_element = ring_element
        for attribute in ring_element.__dict__:
            setattr(self, attribute, ring_element.__dict__[attribute])

    def __repr__(self) -> str:
        return repr(self.ring_element)

    def __hash__(self) -> int:
        return hash(self.ring_element)

    def __mul__(self, other) -> GroupElementAdapter:
        return GroupElementAdapter(self.ring_element * other.ring_element)

    def __eq__(self, other):
        return self.ring_element == other.ring_element

    def __ne__(self, other):
        return self.ring_element != other.ring_element

    def get_order(self):
        return self.ring_element.multiplicative_order()

    def is_identity(self):
        return self.ring_element.is_multiplicative_identity()


class ResidueClass(RingElement):
    """Base class for elements of the ring of integers modulo N."""

    def __init__(self, residue: int, modulus: int):
        self.modulus = modulus
        self.residue = residue % modulus

    def __repr__(self):
        return "[" + str(self.residue) + "]"

    def __hash__(self):
        return hash((self.residue, self.modulus))

    def __eq__(self, other):
        return (self.residue == other.residue) and (self.modulus == other.modulus)

    def __ne__(self, other):
        return not (self.residue == other.residue) or not (
            self.modulus == other.modulus
        )

    def __add__(self, other):
        if self.modulus != other.modulus:
            raise ValueError("the moduli must be equal")
        return ResidueClass((self.residue + other.residue) % self.modulus, self.modulus)

    def __mul__(self, other):
        if self.modulus != other.modulus:
            raise ValueError("the moduli must be equal")
        return ResidueClass((self.residue * other.residue) % self.modulus, self.modulus)

    def __pow__(self, N: int):
        if N > 0:
            return reduce(lambda x, y: x * y, [self] * N)
        if N < 0:
            return reduce(lambda x, y: x * y, [~self] * abs(N))
        return ResidueClass(1, self.modulus)

    def __neg__(self):
        return ResidueClass(-self.residue, self.modulus)

    def __sub__(self, other):
        return self + (-other)

    def __invert__(self):
        if galois.gcd(self.residue, self.modulus) != 1:
            raise ValueError(
                """the residue and modulus must be relatively prime 
                for a multiplicative inverse to exist"""
            )
        # this is Euler's theorem
        exp = galois.euler_phi(self.modulus) - 1
        return ResidueClass(self.residue**exp, self.modulus)

    def additive_order(self):
        return int(self.modulus / galois.gcd(self.residue, self.modulus))

    def multiplicative_order(self):
        prod = self
        i = 1
        while not prod.is_multiplicative_identity():
            prod *= self
            i += 1
        return i

    def is_additive_identity(self) -> bool:
        return self.residue == 0

    def is_multiplicative_identity(self) -> bool:
        return self.residue == 1

    def is_invertible(self) -> bool:
        return galois.gcd(self.residue, self.modulus) == 1
