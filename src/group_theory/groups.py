"""
This module defines the Group, Subgroup, and Coset classes. 
"""

from __future__ import annotations

from random import sample
from itertools import product
from typing import Type

import galois

from group_theory.group_elements import (
    CartesianProductElement,
    GroupElement,
)


class Group:
    """Base class representing a finite group.

    Args:
        elements (list[GroupElement]): a list of elements for the group.
    """

    # pylint: disable=too-many-instance-attributes
    # pylint: disable=too-many-public-methods

    def __init__(self, elements: list[GroupElement]):
        self.elements = elements
        self.order = len(self.elements)
        self.identity = self.get_identity()
        self.canonical_generators = None

        # properties
        self._generators = None
        self._center = None
        self._commutator_subgroup = None
        self._conjugacy_classes = None
        self._generator_representations = None

    def get_orders(self) -> None:
        """Fetches the orders of all elements of the group."""
        for element in self:
            element.get_order()

    def get_identity(self) -> GroupElement:
        """Fetches the identity element of the group

        Returns:
            GroupElement: the identity element of the group.
        """
        for element in self.elements:
            if element.is_identity():
                self.identity = element
                return element
        return None

    def __repr__(self) -> str:
        repr_string = ""
        for group_element in self:
            repr_string += repr(group_element) + "\n"
        return repr_string

    def __mul__(self, other) -> Group:
        """This is the standard cartesian product of groups: for group G and H, G*H contains all
        pairs of elements of the form (g,h)

        Args:
            other (Group): the other Group with which to form the Cartesian product.

        Returns:
            Group: the Cartesian product of the groups self and other.
        """
        product_elements = list(
            CartesianProductElement(g) for g in product(self.elements, other.elements)
        )
        return Group(product_elements)

    def __truediv__(self, other) -> Group:
        """Quotient group

        Args:
            other (Subgroup): a normal subgroup of G(self)

        Raises:
            ValueError: _description_
            TypeError: _description_

        Returns:
            Group: Returns the quotient group G/N for a normal subgroup N(other) of G(self)
        """
        if isinstance(other, Subgroup):
            if other.parent_group == self:
                cosets = {Coset(g, other) for g in self}
                return Group(list(cosets))
            raise ValueError(
                "the subgroup provided is not a subgroup of the provided parent group"
            )
        raise TypeError("you must pass a valid subgroup to quotient by")

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

    def __getitem__(self, key) -> GroupElement:
        return self.elements[key]

    def __iter__(self):
        return iter(self.elements)

    def __len__(self) -> int:
        return len(self.elements)

    def is_trivial(self) -> bool:
        """Determines whether or not the group is trivial, i.e., if G = {1}.

        Returns:
            bool: True if the group is trivial, False otherwise.
        """
        if len(self.elements) == 1 and self.elements[0].is_identity():
            return True
        return False

    def is_abelian(self) -> bool:
        """Determines whether or not the group is abelian, i.e., if x*y=y*x
        for all x,y in G.

        Returns:
            bool: True if the group is abelian, False otherwise.
        """
        return self.commutator_subgroup.is_trivial()

    def is_perfect(self) -> bool:
        """Determines whether or not the group is perfect, i.e., if G = [G,G].

        Returns:
            bool: True if the group is perfect, False otherwise.
        """
        return self.commutator_subgroup == self

    def is_solvable(self) -> bool:
        """Determines whether or not the group is solvable, i.e., if the derived
        series of G(self) ends at the trivial group.

        Returns:
            bool: True if the group is solvable, False otherwise.
        """
        return self.derived_series()[-1].is_trivial()

    def is_nilpotent(self) -> bool:
        """Determines whether or not the group is nilpotent, i.e., if the lower
        central series of G(self) ends at the trivial group.

        Returns:
            bool: True if the group is trivial, False otherwise.
        """
        return self.lower_central_series()[-1].is_trivial()

    def abelianization(self) -> Group:
        """Constructs the abelianization of G, i.e., the quotient group G/[G,G].

        Returns:
            Group: the quotient group G/[G,G].
        """
        return self / self.commutator_subgroup

    def get_random_generators(self) -> list[GroupElement]:
        """Fetches a random generating set for the group.

        Returns:
            list[GroupElement]: a set of generating elements for the group.
        """
        if len(self.elements) == 1:
            return [self[0]]
        generators = []
        elements = [g for g in self if not g.is_identity()]
        generators += [sample(elements, 1)[0]]
        generated_elements = self.subgroup_generated_by(generators).elements
        while len(generated_elements) != self.order:
            elements = [g for g in elements if g not in generated_elements]
            generators += [sample(elements, 1)[0]]
            generated_elements = self.subgroup_generated_by(generators).elements
        return generators

    def get_generators(self) -> list[GroupElement]:
        """If the group has a canonical set of generators, return those. Otherwise,
        find a random set of generators.

        Returns:
            list[GroupElement]: a list of GroupElements which generate the group.
        """
        if self.canonical_generators:
            return self.canonical_generators
        return self.get_random_generators()

    @property
    def generators(self):
        """Fetches a set of generators for the group"""
        if self._generators is None:
            self._generators = self.get_generators()
            return self._generators
        return self._generators

    def get_center(self) -> Subgroup:
        """Computes the center of the group, using the fact that the
        center of a group is equal to the centralizer C_G(G).

        Returns:
            Subgroup: the center of the group.
        """
        return self.centralizer(Subgroup(self.elements, self))

    @property
    def center(self) -> Subgroup:
        """Fetches the center of the group"""
        if self._center is None:
            self._center = self.get_center()
            return self._center
        return self._center

    def get_conjugacy_classes(self):
        """
        each element of the center forms a conjugacy class consisting of just itself,
        so we compute the center first and the other classes second
        """
        conjugacy_classes = set()
        conjugacy_classes.update([frozenset({z}) for z in self.center])
        for h in set(self.elements) - set(self.center.elements):
            conjugacy_classes.add(frozenset({g * h * (~g) for g in self}))
        return conjugacy_classes

    @property
    def conjugacy_classes(self):
        """Fetches the set of conjugacy classes of the group"""
        if self._conjugacy_classes is None:
            self._conjugacy_classes = self.get_conjugacy_classes()
            return self._conjugacy_classes
        return self._conjugacy_classes

    def get_commutator_subgroup(self) -> Subgroup:
        """Computes the commutator subgroup of the group. This is done
        using generators: if G is generated by a set S, then the commutator
        subgroup is the normal closure of the set of commutators of
        elements of S.

        Returns:
            Subgroup: the commutator subgroup of the group.
        """
        if self.generators:
            return self.subgroup_generated_by(
                [
                    g * x * y * (~x) * (~y) * (~g)
                    for x in self.generators
                    for y in self.generators
                    for g in self.elements
                ]
            )
        return self.subgroup_generated_by(
            [x * y * (~x) * (~y) for x in self.elements for y in self.elements]
        )

    @property
    def commutator_subgroup(self):
        """Fetches the commutator subgroup"""
        if self._commutator_subgroup is None:
            self._commutator_subgroup = self.get_commutator_subgroup()
            return self._commutator_subgroup
        return self._commutator_subgroup

    def centralizer(self, subgroup: Subgroup) -> Subgroup:
        """Computes the centralizer of a subgroup.

        Args:
            H (Subgroup): A subgroup of of the group (self)

        Raises:
            ValueError: if the parent group of the subgroup is not (self), then
            this error is raised.

        Returns:
            Subgroup: the centralizer subgroup of the subgroup passed.
        """
        if self != subgroup.parent_group:
            raise ValueError("subgroup provided is not a subgroup of this group")
        centralizer_elements = [self.identity]
        number_of_generators = len(subgroup.generators)
        for z in [x for x in self if x != self.identity]:
            counter = 0
            for h in subgroup.generators:
                if z * h == h * z:
                    counter += 1
                    continue
                else:
                    break
            if counter == number_of_generators:
                centralizer_elements += [z]
        return Subgroup(centralizer_elements, self)

    def normalizer(self, subgroup: Subgroup) -> Subgroup:
        """Computes the normalizer of a subgroup.

        Args:
            subgroup (Subgroup): the subgroup whose normalizer we want to compute.

        Raises:
            ValueError: if the parent group of the subgroup is not (self), then
            this error is raised.

        Returns:
            Subgroup: the normalizer of the subgroup provided.
        """
        if self != subgroup.parent_group:
            raise ValueError("subgroup provided is not a subgroup of this group")
        normalizer_elements = [self.identity]
        for z in [x for x in self if x != self.identity]:
            if subgroup.left_coset(z) == subgroup.right_coset(z):
                normalizer_elements += [z]
        return Subgroup(normalizer_elements, self)

    def subgroup_generated_by(self, generators: list[GroupElement]) -> Subgroup:
        """Computes the subgroup generated by a list of elements of the group.

        Args:
            generators (list[GroupElement]): the list of generators for the subgroup.

        Raises:
            TypeError: this error is raised if a list is not passed.

        Returns:
            Subgroup: the subgroup of (self) generated by the elements provided.

        Note: this algorithm is the "black-box" algo found here:
        https://groupprops.subwiki.org/w/index.php?title=Black-box_group_algorithm_for_finding_the_subgroup_generated_by_a_subset
        """
        if not isinstance(generators, list):
            raise TypeError("the generators must be a list of group elements")
        generated_elements = {self.identity}
        length_i_words = {self.identity}

        while length_i_words:
            products = {f * s for f in length_i_words for s in generators}
            products = products - generated_elements
            length_i_words = products
            generated_elements = generated_elements.union(length_i_words)

        return Subgroup(list(generated_elements), self)

    def get_generator_representations(self) -> dict:
        """Fetches the generator representations for each element of the group.
        For example, if G has generators g_1 and g_2, and g = (g_1^3)*(g_2^2)*(g_1), then
        the list associated with g is [g_1, g_1, g_1, g_2, g_2, g_1].

        Returns:
            dict: a dictionary whose keys are elements of the group, with values
            lists whose elements are generators of G, with product equal to g (the key).
        """
        generated_elements = {self.identity}
        length_i_words = {self.identity}
        generator_dict = {self.identity: [self.identity]}
        for g in self.generators:
            generator_dict[g] = [g]

        while length_i_words:
            products = set()
            for word in length_i_words:
                for generator in self.generators:
                    x = word * generator
                    products.add(x)
                    if x not in generator_dict:
                        # note: since F={e} to start, f is always in generator_dict.keys()
                        generator_dict[x] = generator_dict[word] + [generator]
            products = products - generated_elements
            length_i_words = products
            generated_elements = generated_elements.union(length_i_words)

        return generator_dict

    @property
    def generator_representations(self):
        """Fetches generator representations for each element"""
        if self._generator_representations is None:
            self._generator_representations = self.get_generator_representations()
            return self._generator_representations
        return self._generator_representations

    def derived_series(self) -> list[Subgroup]:
        """Computes the derived series of the group.

        Returns:
            list: a list of Subgroups which make up the derived series of the group.
        """
        last_subgroup = Subgroup(self.elements, self)
        derived_series_list = [last_subgroup]
        next_subgroup = last_subgroup.commutator_subgroup
        while next_subgroup != last_subgroup:
            derived_series_list.append(next_subgroup)
            last_subgroup = next_subgroup
            next_subgroup = last_subgroup.commutator_subgroup
        return derived_series_list

    def lower_central_series(self) -> list[Subgroup]:
        """Computes the lower central series of the group.

        Returns:
            list: a list of Subgroups which make up the lower central series of the group.
        """
        last_subgroup = Subgroup(self.elements, self)
        lower_central_series_list = [last_subgroup]
        next_subgroup = last_subgroup.commutator_subgroup
        while next_subgroup != last_subgroup:
            lower_central_series_list.append(next_subgroup)
            last_subgroup = next_subgroup
            next_subgroup = self.subgroup_generated_by(
                [x * y * (~x) * (~y) for x in last_subgroup for y in self]
            )
        return lower_central_series_list

    def upper_central_series(self) -> list[Subgroup]:
        """Computes the upper central series of the group.

        Returns:
            list: a list of Subgroups which make up the upper central series of the group.
        """
        last_center = Subgroup([self.identity], self)
        upper_central_series_list = [last_center]
        next_center = self.center
        while last_center != next_center:
            upper_central_series_list.append(next_center)
            last_center = next_center
            next_center = Subgroup(
                [
                    x
                    for x in self
                    if {x * y * (~x) * (~y) for y in self}.issubset(last_center)
                ],
                self,
            )
        return upper_central_series_list


class LinearGroup(Group):
    """
    Base class for groups of invertible matrices over finite fields.

    In the future, this class will be used for group representations.
    """

    def __init__(
        self, elements: list[GroupElement], ground_field: Type[galois.FieldArray]
    ):
        super().__init__(elements)
        self.ground_field = ground_field


class Subgroup(Group):
    """Base class for subgroups of finite groups.

    Args:
        - parent_group (Group): the group which the subgroup is a subset of.
        - elements (list): the list of elements of the subgroup.
    """

    def __init__(self, elements: list[GroupElement], parent_group: Group):
        super().__init__(elements)
        self.parent_group = parent_group
        self.canonical_generators = None
        self.index = int(self.parent_group.order / self.order)

        # properties
        self._is_normal = None

    def validate_inclusion(self) -> bool:
        """Checks whether the subgroup is a subset of the parent group.

        Returns:
            bool: returns True if the subgroup is a subset of the parent group,
            False otherwise.
        """
        return set(self.elements).issubset(self.parent_group.elements)

    # IMPORTANT - until this is called and returns true,
    # the Subgroup instance may not actually be a subgroup
    def validate_subgroup(self) -> bool:
        """Cehcks whether or not the subgroup satisfies the group axioms.

        Returns:
            bool: returns True if the subgroup satisfies the group axioms, returns
            False otherwise.
        """
        for g in self:
            for h in self:
                if g * (~h) not in self:
                    return False
        return True

    def __repr__(self):
        return "Subgroup(" + str(self.generators) + ")"

    def __eq__(self, other):
        return set(self.elements) == set(other.elements)

    def __hash__(self):
        return hash((self.elements, self.parent_group))

    def __and__(self, other):
        if self.parent_group != other.parent_group:
            raise ValueError("the subgroups must both be a subset of the same group")
        return Subgroup(
            [x for x in self.elements if x in other.elements], self.parent_group
        )

    def __matmul__(self, other):
        if self.parent_group != other.parent_group:
            raise ValueError("the subgroups must both be a subset of the same group")
        right_product = {h * k for h in self for k in other}
        left_product = {k * h for h in self for k in other}
        if right_product == left_product:
            return Subgroup(list(right_product), self.parent_group)
        else:
            return right_product

    def left_coset(self, left_factor: GroupElement) -> set:
        """Computes the left coset of the subgroup for the given left factor.

        Args:
            left_factor (GroupElement): the group element which will act on all
            subgroup elements on the left

        Raises:
            ValueError: this error is raised if the left_factor provided is not
            an element of the parent group.

        Returns:
            set: the left coset of the subgroup for the given left factor.
        """
        if left_factor not in self.parent_group:
            raise ValueError("group element must be a member of the parent group")
        return {left_factor * x for x in self.elements}

    def right_coset(self, right_factor: GroupElement) -> set:
        """Computes the right coset of the subgroup for the given left factor.

        Args:
            left_factor (GroupElement): the group element which will act on all
            subgroup elements on the right

        Raises:
            ValueError: this error is raised if the right_factor provided is not
            an element of the parent group.

        Returns:
            set: the right coset of the subgroup for the given right factor.
        """
        if right_factor not in self.parent_group:
            raise ValueError("group element must be a member of the parent group")
        return {x * right_factor for x in self.elements}

    def conjugate_subgroup(self, g: GroupElement) -> Subgroup:
        """Computes the conjugate subgroup for a given element of the parent group.
        If H is the subgroup in question, then this method returns the subgroup gHg^{-1}.

        Args:
            g (GroupElement): the element of the parent group to conjugate everything by.

        Raises:
            ValueError: this error is raised if the GroupElement provided is not
            an element of the parent group.

        Returns:
            Subgroup: the conjugate subgroup g(self)g^{-1}.
        """
        if g not in self.parent_group:
            raise ValueError("group element must be a member of the parent group")
        return Subgroup([g * x * (~g) for x in self], self.parent_group)

    def check_normality(self) -> bool:
        """Checks whether or not the subgroup is a normal subgroup of the parent group.

        Returns:
            bool: returns True if the subgroups is normal, returns False otherwise.
        """
        if self.index == 2:  # index 2 subgroups are always normal
            return True
        for g in self.parent_group.generators:
            if self.left_coset(g) != self.right_coset(g):
                return False
        return True

    @property
    def is_normal(self) -> bool:
        """Fetches the normality boolean value."""
        if self._is_normal is None:
            self._is_normal = self.check_normality()
            return self._is_normal
        return self._is_normal


class Coset(GroupElement):
    """
    this class is used to represent elements gH of quotient groups G/H
    for elements g and normal subgroups H
    """

    def __init__(self, g: GroupElement, subgroup: Subgroup) -> None:
        self.g, self.subgroup = self.validate_coset(g, subgroup)
        self.elements = [self.g * x for x in self.subgroup]

    def __repr__(self):
        return "Coset(" + str(self.g) + ")"

    def __eq__(self, other):
        return self.g * (~other.g) in self.subgroup

    def __ne__(self, other):
        return self.g * (~other.g) not in self.subgroup

    def __hash__(self):
        return hash(frozenset(self.elements))

    def __getitem__(self, key):
        return self.elements[key]

    def __iter__(self):
        return iter(self.elements)

    def __len__(self):
        return len(self.elements)

    def __mul__(self, other):
        if self.subgroup != other.subgroup:
            raise ValueError("the subgroups of each coset must match")
        return Coset(self.g * other.g, self.subgroup)

    def __invert__(self):
        return Coset(~self.g, self.subgroup)

    @staticmethod
    def validate_coset(g: GroupElement, subgroup: Subgroup):
        """Checks to make sure that the left factor is an element of the parent group.

        Args:
            g (GroupElement): the left factor
            subgroup (Subgroup): the subgroup to form the coset out of

        Returns:
            g, subgroup (GroupElement, Subgroup): returns the GroupElement and Subgroup if
            the check passes.
        """
        assert (
            g in subgroup.parent_group
        ), "the group element must lie in the parent group"
        return g, subgroup

    def get_order(self):
        return NotImplemented

    def is_identity(self):
        return self.g in self.subgroup
