"""This module defines the GroupAction class."""

from typing import (
    Callable,
    Iterable,
    Tuple,
    Type,
)

from group_theory.group_elements import (
    GroupElement,
)
from group_theory.groups import (
    Group,
    Subgroup,
)


class GroupAction:
    """Base class for an action of a finite group G on a finite set X."""

    def __init__(
        self, G: Group, X: Iterable, action: Callable[[Tuple[GroupElement, Type]], Type]
    ) -> None:
        self.G = G
        self.X = X
        self.action = action

        # properties
        self._orbits = None
        self._fixed_points = None
        self._is_transitive = None
        self._is_faithful = None
        self._is_free = None

    def validate_action(self) -> bool:
        """Verifies that the group action axioms are satisfied.

        Returns:
            bool: returns True if the axioms for a group action are satisfied,
            returns False otherwise.
        """
        identity_condition = all((self.action(self.G.identity, x) == x for x in self.X))
        compatibility_condition = []
        for g in self.G:
            for h in self.G:
                compatibility_condition.append(
                    all(
                        (
                            self.action(g, self.action(h, x)) == self.action(g * h, x)
                            for x in self.X
                        )
                    )
                )
        return identity_condition and all(compatibility_condition)

    def stabilizer(self, x: Type) -> Subgroup:
        """Computes the stabilizer of an element.

        Args:
            x (Type): the element whose stabilizer we want to compute

        Returns:
            Subgroup: the stabilizer subgroup of the element provided.
        """
        return Subgroup([g for g in self.G if self.action(g, x) == x], self.G)

    def get_orbits(self) -> set:
        """Computes the set of orbits of the action.

        Returns:
            set: the set of orbits of the action.
        """
        orbit_set = set()
        for x in self.X:
            orbit_set.add(frozenset({self.action(g, x) for g in self.G}))
        return orbit_set

    @property
    def orbits(self) -> set:
        """Fetches the orbits of the action"""
        if self._orbits is None:
            self._orbits = self.get_orbits()
            return self._orbits
        return self._orbits

    def get_fixed_points(self) -> set:
        """Computes the set of fixed points of the group action."""
        return {
            x
            for x in self.X
            if all((self.action(g, x) == x for g in self.G.generators))
        }

    @property
    def fixed_points(self) -> set:
        """Fetches the fixed points of the action"""
        if self._fixed_points is None:
            self._fixed_points = self.get_fixed_points()
            return self._fixed_points
        return self._fixed_points

    def check_transitivity(self) -> bool:
        """Determines whether or not the group action is transitive."""
        return len(self.orbits) == 1

    @property
    def is_transitive(self) -> bool:
        """Fetches the is_transitive property"""
        if self._is_transitive is None:
            self._is_transitive = self.check_transitivity()
            return self._is_transitive
        return self._is_transitive

    def check_faithfulness(self) -> bool:
        """Determines whether or not the group action is faithful."""
        return (
            len([g for g in self.G if all((self.action(g, x) == x for x in self.X))])
            == 1
        )

    @property
    def is_faithful(self):
        """Fetches the is_faithful property"""
        if self._is_faithful is None:
            self._is_faithful = self.check_faithfulness()
            return self._is_faithful
        return self._is_faithful

    def check_if_free(self):
        """Determines whether or not the group action is free."""
        return (
            len([g for g in self.G if any((self.action(g, x) == x for x in self.X))])
            == 1
        )

    @property
    def is_free(self):
        """Fetches the is_free property"""
        if self._is_free is None:
            self._is_free = self.check_if_free()
            return self._is_free
        return self._is_free
