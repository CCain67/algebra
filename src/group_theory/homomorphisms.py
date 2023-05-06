"""
This module defines the Homomorphism, GroupHom, and Automorphism classes,
and related functions.
"""

from __future__ import annotations
from typing import Callable, Tuple, Union
from functools import reduce
from itertools import product
from copy import copy

from group_theory.group_elements import (
    CartesianProductElement,
    GroupElement,
)
from group_theory.groups import (
    Group,
    Subgroup,
)


def _homomorphism_factory(in_out_dict: dict) -> Callable:
    def homomorphism(g: GroupElement) -> GroupElement:
        return in_out_dict[g]

    return homomorphism


class Homomorphism:
    """Base class for homomorphisms between fininte groups.

    Args:
        - domain (Group): the group of inputs of the homomorphism
        - morphism (Callable): the function/map defining the homomorphism
        - codomain (Group): the group of outputs of the homomorphism.
    """

    def __init__(
        self,
        domain: Group,
        morphism: Callable[[Tuple[GroupElement, ...]], GroupElement],
        codomain: Group,
    ) -> None:
        self.domain = domain
        self.morphism = morphism
        self.codomain = codomain
        self.in_out_dict = None

        # properties
        self._graph = None
        self._image = None
        self._kernel = None
        self._is_iso = None

    def __repr__(self) -> str:
        repr_string = ""
        for g in self.domain.generators:
            repr_string += str(g) + " -> " + str(self.morphism(g)) + "\n"
        return repr_string

    def __mul__(self, other) -> Homomorphism:
        if not isinstance(other, Homomorphism):
            raise TypeError("type must be Homomorphism")
        if self.domain != other.codomain:
            raise ValueError("non-composable homomorphisms")
        self._morphism_to_in_out_dict()
        in_out_dict = {g: self(other(g)) for g in other.domain.generators}
        return Homomorphism.from_action_on_generators(
            other.domain, in_out_dict, self.codomain
        )

    def __pow__(self, power: int) -> Homomorphism:
        if self.domain != self.codomain:
            raise ValueError(
                "Composing a Homomorphism with itself requires the domain and codomain to be equal"
            )
        if power > 0:
            return reduce(lambda x, y: x * y, [self] * power)
        if power < 0:
            return reduce(lambda x, y: x * y, [~self] * abs(power))
        return Homomorphism.from_dict(
            self.domain, {g: g for g in self.domain}, self.codomain
        )

    def __invert__(self) -> Homomorphism:
        if not self.is_iso:
            raise ValueError("only isomorphisms may be inverted")
        return Homomorphism.from_dict(
            self.codomain, {self(g): g for g in self.domain}, self.codomain
        )

    def __eq__(self, other) -> bool:
        if self.domain != other.domain:
            raise ValueError("the domains of the homomorphisms are not equal")
        if self.codomain != other.codomain:
            raise ValueError("the codomains of the homomorphisms are not equal")
        return all(
            (self.morphism(g) == other.morphism(g) for g in self.domain.generators)
        )

    def __ne__(self, other) -> bool:
        if self.domain != other.domain:
            raise ValueError("the domains of the homomorphisms are not equal")
        if self.codomain != other.codomain:
            raise ValueError("the codomains of the homomorphisms are not equal")
        return any(
            (self.morphism(g) != other.morphism(g) for g in self.domain.generators)
        )

    def __hash__(self) -> int:
        in_out_generators = {g: self(g) for g in self.domain.generators}
        return hash((self.domain, frozenset(in_out_generators.items()), self.codomain))

    def __call__(self, g) -> GroupElement:
        try:
            return self.morphism(g)
        except KeyError:
            self._fetch_remaining_images()
            return self.morphism(g)

    def is_identity(self) -> bool:
        """Determines whether or not the homomorphism is the identity.

        Returns:
            bool: returns True if the homomorphism is the identity map, returns False otherwise.
        """
        return all(
            (
                self.domain == self.codomain,
                *[self.morphism(g) == g for g in self.domain.generators],
            )
        )

    def is_trivial(self) -> bool:
        """Determines whether or not the homomorphism is trivial.

        Returns:
            bool: returns True if the homomorphism is trivial, returns False otherwise.
        """
        return all(
            (self.morphism(g) == self.codomain.identity for g in self.domain.generators)
        )

    def validate_homomorphism(self) -> bool:
        """Verifies that the homomorphism is actually a homomorphism.

        Returns:
            bool: returns True if the homomorphism is actually a homomorphism,
            returns False otherwise.

        Explanation:
            this comes from the paper 'Computing with Group Homomorphisms' by Leedham-Green et. al.

            The main idea is that a map f:G -> H is a homomorphism if and only if the subgroup
            generated by (g,f(g)) (see the method get_graph below) has trivial intersection
            with {1}xH in GxH, if and only if this generated subgroup has the same order as G.
        """
        return self.graph.order == self.domain.order

    def _fetch_remaining_images(self) -> None:
        """In some cases, e.g. when defining Hom-sets, homomorphisms are defined only on
        generators for efficiency, and so the rest of the images need to be determined.
        """
        generator_in_out_dict = {g: self.morphism(g) for g in self.domain.generators}

        generator_in_out_dict[self.domain.identity] = self.codomain.identity
        generator_representations = self.domain.generator_representations
        in_out_dict = {
            g: reduce(
                lambda x, y: x * y,
                [generator_in_out_dict[x] for x in generator_representations[g]],
            )
            for g in self.domain
        }
        self.morphism = _homomorphism_factory(in_out_dict=in_out_dict)

    def _morphism_to_in_out_dict(self) -> None:
        try:
            in_out_dict = {g: self.morphism(g) for g in self.domain}
        except KeyError:
            self._fetch_remaining_images()
            in_out_dict = {g: self.morphism(g) for g in self.domain}
        self.in_out_dict = in_out_dict
        self.morphism = _homomorphism_factory(in_out_dict=in_out_dict)

    @classmethod
    def from_action_on_generators(
        cls, domain: Group, generator_in_out_dict: dict, codomain: Group
    ) -> Homomorphism:
        """Defines a homomorphism via its action on generators.
        TODO: add validation

        Args:
            domain (Group): the group of inputs
            generator_in_out_dict (dict): a dictionary of the form {generator:f(generator)}
            codomain (Group): the group of outputs

        Raises:
            ValueError: this error is raised if the generators do not form a subset of the
            domain provided.

        Returns:
            Homomorphism: a homomorphism defined on the generators provided.
        """
        if not set(domain.generators).issubset(set(generator_in_out_dict.keys())):
            raise ValueError(
                "the generators of the domain must be present in the dictionary keys"
            )

        generator_in_out_dict[domain.identity] = codomain.identity
        generator_representations = domain.generator_representations
        in_out_dict = {
            g: reduce(
                lambda x, y: x * y,
                [generator_in_out_dict[x] for x in generator_representations[g]],
            )
            for g in domain
        }

        return cls.from_dict(domain=domain, in_out_dict=in_out_dict, codomain=codomain)

    @classmethod
    def from_dict(
        cls, domain: Group, in_out_dict: dict, codomain: Group
    ) -> Homomorphism:
        """Creates a homomorphism from a input-output dictionary

        Args:
            domain (Group): group of inputs
            in_out_dict (dict): a dictionary of the form {input:output}
            codomain (Group): group of outputs

        Returns:
            Homomorphism: a homomorphism F which satisfies F(input) = in_out_dict[input].
        """
        pre_homomorphism = _homomorphism_factory(in_out_dict=in_out_dict)
        homomorphism = cls(domain, pre_homomorphism, codomain)
        homomorphism.in_out_dict = in_out_dict
        return homomorphism

    def get_graph(self) -> Subgroup:
        """
        This is the subgroup of the cartesian product GxH generated by (g,f(g)) for a map f:G -> H,
        mentioned in the docstring above.
        """
        augmented_generators = self.domain.generators + [
            a * b for a in self.domain.generators for b in self.domain.generators
        ]
        cartesian_product = self.domain * self.codomain
        element_image_pairs = [
            CartesianProductElement((x, self.morphism(x))) for x in augmented_generators
        ]
        graph = cartesian_product.subgroup_generated_by(element_image_pairs)
        return graph

    @property
    def graph(self) -> Subgroup:
        """Fetches the graph subgroup of the cartesian product domain*codomain."""
        if self._graph is None:
            self._graph = self.get_graph()
            return self._graph
        return self._graph

    def get_image(self) -> Subgroup:
        """Computes the image of the homomorphism.

        Returns:
            Subgroup: the image of the homomorphism as a subgroup of the codomain.
        """
        image_generators = [self.morphism(g) for g in self.domain.generators]
        return self.codomain.subgroup_generated_by(image_generators)

    @property
    def image(self) -> Subgroup:
        """Fetches the image of the homomorphism."""
        if self._image is None:
            self._image = self.get_image()
            return self._image
        return self._image

    def get_kernel(self) -> Subgroup:
        """Computes the kernel of the homomorphism.

        Returns:
            Subgroup: the kernel of the homomorphism as a subgroup of the domain.
        """
        if self.in_out_dict:
            if len(self.in_out_dict) != len(self.domain):
                self._fetch_remaining_images()
                return Subgroup(
                    [
                        g
                        for g in self.domain
                        if self.in_out_dict[g] == self.codomain.identity
                    ],
                    self.domain,
                )
            return Subgroup(
                [
                    g
                    for g in self.domain
                    if self.in_out_dict[g] == self.codomain.identity
                ],
                self.domain,
            )
        graph = self.graph
        cartesian_product = self.domain * self.codomain
        domain_times_identity = (
            Subgroup(  # this is the subgroup domain x {1} of domain x codomain
                [p for p in cartesian_product if p[1].is_identity()], cartesian_product
            )
        )
        return Subgroup([p[0] for p in (graph & domain_times_identity)], self.domain)

    @property
    def kernel(self) -> Subgroup:
        """Fetches the kernel of the homomorphism"""
        if self._kernel is None:
            self._kernel = self.get_kernel()
        return self._kernel

    def check_iso(self) -> bool:
        """Determines if the homomorphism is an isomorphism.

        Returns:
            bool: returns True if the homomorphism is an isomorphism, returns False otherwise.

        Note:
            This utilizes the fact that for two finite sets, if an injective map exists between
            them, and the two sets are of equal order, then the map is a bijection.
        """
        if (
            len(self.kernel) == 1
            and self.domain.order == self.codomain.order
            and self.validate_homomorphism()
        ):
            return True
        return False

    @property
    def is_iso(self) -> bool:
        """Fetches the is_iso property."""
        if self._is_iso is None:
            self._is_iso = self.check_iso()
        return self._is_iso


class GroupHom:
    """
    Class representing the hom-set Hom(G,H) for two finite groups G and H.
    """

    def __init__(self, domain: Group, codomain: Group) -> None:
        self.domain = domain
        self.codomain = codomain
        self.homomorphisms = self.get_all_homomorphisms()

    def __repr__(self) -> str:
        return str(self.homomorphisms)

    def __getitem__(self, key) -> Homomorphism:
        return self.homomorphisms[key]

    def __iter__(self):
        return iter(self.homomorphisms)

    def __hash__(self) -> int:
        return hash(frozenset(self.homomorphisms))

    def __len__(self) -> int:
        return len(self.homomorphisms)

    def __eq__(self, other) -> bool:
        return set(self.homomorphisms) == set(other.homomorphisms)

    def __ne__(self, other) -> bool:
        return set(self.homomorphisms) != set(other.homomorphisms)

    def get_all_homomorphisms(self) -> list:
        """
        fetching all possible homomorphisms works as follows:

        1. for each generator g of G, we find all elements h in H such that
        order(h) divides order(g) - these are the candidates for the image of g

        2. if f:G -> H is a homomorphism, we must have f(xy)=f(x)f(y), so we
        check that order(f(x)f(y)) divides order of xy. if not, we move on to the
        next set of potential images of generators.

        3. once we are past this check, we call validate_homomorphism. if this passes,
        we append the homomorphism to the list.
        """
        possible_images = {
            g: [h for h in self.codomain if g.order % h.order == 0]
            for g in self.domain.generators
        }

        possible_homomorphisms = product(
            *[possible_images[g] for g in possible_images.keys()]
        )

        pre_hom_set = [
            {g: hom_base[i] for i, g in enumerate(possible_images.keys())}
            for hom_base in possible_homomorphisms
        ]

        homset = []
        for in_out_dict in pre_hom_set:
            go_to_next_in_out = False
            for a, b in product(self.domain.generators, self.domain.generators):
                if (a * b).order % (in_out_dict[a] * in_out_dict[b]).order == 0:
                    in_out_dict[a * b] = in_out_dict[a] * in_out_dict[b]
                else:
                    go_to_next_in_out = True
                    break
            if go_to_next_in_out:
                continue
            homomorphism = Homomorphism(
                domain=self.domain,
                morphism=_homomorphism_factory(in_out_dict),
                codomain=self.codomain,
            )
            if homomorphism.validate_homomorphism():
                homset.append(homomorphism)
        return homset


class Automorphism(Homomorphism, GroupElement):
    """Class representing an automorphism of a finite group."""

    def __init__(
        self, group: Group, morphism: Callable[[Tuple[GroupElement, ...]], GroupElement]
    ) -> None:
        super().__init__(group, morphism, group)
        self.group = group

        # properties
        self._order = None

    def __mul__(self, other) -> Union[Homomorphism, Automorphism]:
        if not isinstance(other, Homomorphism):
            raise TypeError("type must be Homomorphism")
        if self.domain != other.codomain:
            raise ValueError("non-composable homomorphisms")
        if isinstance(other, Automorphism):
            self._morphism_to_in_out_dict()
            in_out_dict = {g: self(other(g)) for g in other.domain.generators}
            return Automorphism.from_action_on_generators(
                G=self.group, generator_in_out_dict=in_out_dict
            )
        self._morphism_to_in_out_dict()
        in_out_dict = {g: self(other(g)) for g in other.domain.generators}
        return Homomorphism.from_action_on_generators(
            other.domain, in_out_dict, self.codomain
        )

    def __pow__(self, N: int) -> Automorphism:
        if N > 0:
            return reduce(lambda x, y: x * y, [self] * N)
        if N < 0:
            return reduce(lambda x, y: x * y, [~self] * abs(N))
        return Automorphism.from_dict(self.domain, {g: g for g in self.domain})

    def __invert__(self) -> Automorphism:
        return Automorphism.from_dict(self.domain, {self(g): g for g in self.domain})

    def get_order(self) -> int:
        prod = self
        i = 1
        while not prod.is_identity():
            prod *= self
            i += 1
        return i

    @property
    def order(self) -> int:
        """Fetches the order of the automorphism."""
        if self._order is None:
            self._order = self.get_order()
            return self._order
        return self._order

    @classmethod
    def from_action_on_generators(
        cls, G: Group, generator_in_out_dict: dict
    ) -> Automorphism:
        """Creates an automorphism based on inputs and outputs of group generators.

        Args:
            G (Group): input group
            generator_in_out_dict (dict): a dictionary of the form {input:output}

        Raises:
            ValueError: this error is raised if the generators do not form a subset
            of the domain.

        Returns:
            Automorphism: an automorphism defined from the action on the generators
        """
        if not set(G.generators).issubset(set(generator_in_out_dict.keys())):
            raise ValueError(
                "the generators of the domain must be present in the dictionary keys"
            )

        generator_in_out_dict[G.identity] = G.identity
        generator_representations = G.generator_representations
        in_out_dict = {
            g: reduce(
                lambda x, y: x * y,
                [generator_in_out_dict[x] for x in generator_representations[g]],
            )
            for g in G
        }

        return cls.from_dict(G=G, in_out_dict=in_out_dict)

    @classmethod
    def from_dict(cls, G: Group, in_out_dict: dict) -> Automorphism:
        """Creates an automorphism from a input-output dictionary

        Args:
            input_group (Group): group of inputs and outputs
            in_out_dict (dict): a dictionary of the form {input:output}

        Returns:
            Automorphism: an automorphism F which satisfies F(input) = in_out_dict[input].
        """
        pre_automorphism = _homomorphism_factory(in_out_dict=in_out_dict)
        automorphism = cls(G, pre_automorphism)
        automorphism.in_out_dict = in_out_dict
        return automorphism


def _aut_get_identity(group: Group) -> Automorphism:
    identity_dict = {g: g for g in group}
    eye = Automorphism.from_dict(G=group, in_out_dict=identity_dict)
    eye.in_out_dict = identity_dict
    return eye


def _inner_automorphism_factory(g: GroupElement) -> Callable:
    def inner_auto(x: GroupElement) -> GroupElement:
        return g * x * (~g)

    return inner_auto


def Inn(group: Group) -> Group:
    """Creates the group Inn(G) of inner automorphisms of the finite group G."""
    inner_automorphisms = list(
        {Automorphism(group, _inner_automorphism_factory(g)) for g in group}
    )
    inner_automorphism_group = Group(inner_automorphisms)
    inner_automorphism_group.identity = _aut_get_identity(group)
    return inner_automorphism_group


def _check_order_condition(group: Group, in_out_dict: dict) -> bool:
    for a, b in product(group.generators, group.generators):
        order_condition = (a * b).order == (in_out_dict[a] * in_out_dict[b]).order
        if not order_condition:
            return False
    return True


def _check_preserves_subgroup(
    automorphism: Automorphism, subgroup: Subgroup, needs_verification: bool = False
) -> bool:
    if needs_verification:
        return all((automorphism(h) in subgroup for h in subgroup.generators))
    return True


def _check_class_preserving(
    automorphism: Automorphism, needs_verification: bool = False
) -> bool:
    if needs_verification:
        conjugacy_classes_preserved = []
        for conjugacy_class in automorphism.group.conjugacy_classes:
            class_was_preserved = all(
                (automorphism(g) in conjugacy_class for g in conjugacy_class)
            )
        conjugacy_classes_preserved.append(class_was_preserved)

        return all(conjugacy_classes_preserved)
    return True


def Aut(
    group: Group, relative_subgroup: Subgroup = None, class_preserving: bool = False
) -> Group:
    """
    Creates the group Aut(G) of automorphisms of the finite group G.

    Args:
        - G: Group which we are forming the automorphism group of.
        - relative_subgroup: if provided a relative_subgroup H, this function
        computes the subgroup of Aut(G) consisting of automorphisms f:G -> G
        such that f(H) = H.
        - class_preserving: If True, this function returns the subgroup of all
        automorphisms of G which preserve conjugacy classes.

    Fetching all possible automorphisms works just as above in GroupHomSet,
    but with some slight modifications:

    0. First, construct all of the inner automorphisms, since it is not efficient
    to have to discover them all over again.

    1. For each generator g of G, we find all elements h in H such that
    order(h) = order(g) - these are the candidates for the image of g

    2. If f:G -> H is an automorphism, we must have f(xy)=f(x)f(y)
    and order(f(x)f(y)) = order(xy). If not, we move on to the
    next set of potential images of generators.

    3. Once we are past this check, we call validate_homomorphism and check if the
    kernel is trivial. If this passes, we append the automorphism to the list.
    """
    automorphisms = []
    verify_preserves_subgroup = relative_subgroup is not None

    possible_generator_images = {
        g: [h for h in group if g.order == h.order] for g in group.generators
    }
    generator_image_choices = product(
        *[possible_generator_images[g] for g in possible_generator_images.keys()]
    )
    potential_automorphisms = [
        {g: image_list[i] for i, g in enumerate(possible_generator_images.keys())}
        for image_list in generator_image_choices
    ]

    inner_automorphisms = list(
        {
            Automorphism(group, _inner_automorphism_factory(g))
            for g in [group.identity] + [x for x in group if x not in group.center]
        }
    )
    if verify_preserves_subgroup:
        inner_automorphisms = [
            inn
            for inn in inner_automorphisms
            if all((inn(h) in relative_subgroup for h in relative_subgroup.generators))
        ]
    automorphisms += inner_automorphisms
    inner_auto_generator_images = [
        {g: inn(g) for g in group.generators} for inn in inner_automorphisms
    ]

    potential_automorphisms = [
        in_out_dict
        for in_out_dict in potential_automorphisms
        if in_out_dict not in inner_auto_generator_images
    ]
    potential_automorphisms = [
        in_out_dict
        for in_out_dict in potential_automorphisms
        if _check_order_condition(group, in_out_dict)
    ]
    updated_potential_automorphisms = copy(potential_automorphisms)

    for in_out_dict in potential_automorphisms:
        if in_out_dict in updated_potential_automorphisms:
            automorphism = Automorphism.from_action_on_generators(group, in_out_dict)

            homomorphism_check = automorphism.validate_homomorphism()
            kernel_check = len(automorphism.kernel) == 1
            preserves_conjugacy_classes = _check_class_preserving(
                automorphism, class_preserving
            )
            fixes_subgroup = _check_preserves_subgroup(
                automorphism, relative_subgroup, verify_preserves_subgroup
            )

            if all(
                [
                    homomorphism_check,
                    kernel_check,
                    preserves_conjugacy_classes,
                    fixes_subgroup,
                ]
            ):
                automorphisms = list(
                    {A * automorphism for A in automorphisms}
                    .union({automorphism * A for A in automorphisms})
                    .union(set(automorphisms))
                )
                generator_images_to_remove = [
                    {g: A(g) for g in group.generators} for A in automorphisms
                ]
                updated_potential_automorphisms = [
                    in_out_dict
                    for in_out_dict in updated_potential_automorphisms
                    if in_out_dict not in generator_images_to_remove
                ]
            else:
                generator_images_to_remove = [
                    {g: A(in_out_dict[g]) for g in in_out_dict.keys()}
                    for A in automorphisms
                ]
                updated_potential_automorphisms = [
                    in_out_dict
                    for in_out_dict in updated_potential_automorphisms
                    if in_out_dict not in generator_images_to_remove
                ]
        else:
            continue

    automorphism_group = Group(automorphisms)
    automorphism_group.identity = _aut_get_identity(group)
    return automorphism_group


def Out(group: Group) -> Group:
    """
    Creates the group Out(G) of outer automorphisms of the finite group G.

    Out(G) is defined as the quoteint group Aut(G)/Inn(G).
    """
    automorphism_group = Aut(group)
    inner_automorphism_group = Inn(group)
    return automorphism_group / Subgroup(
        inner_automorphism_group.elements, automorphism_group
    )
