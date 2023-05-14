"""
This module contains functions for defining various automorphism groups: 
Inn(G), Aut(G), and Out(G)
"""

from __future__ import annotations
from typing import Callable
from itertools import product
from copy import copy

from group_theory.group_elements import GroupElement
from group_theory.groups import (
    Group,
    Subgroup,
)
from group_theory.homomorphisms import Automorphism


def _aut_get_identity(group: Group) -> Automorphism:
    identity_dict = {g: g for g in group}
    eye = Automorphism.from_dict(G=group, in_out_dict=identity_dict)
    eye.in_out_dict = identity_dict
    return eye


def _inner_automorphism_factory(g: GroupElement) -> Callable:
    def inner_auto(x: GroupElement) -> GroupElement:
        return g * x * (~g)

    return inner_auto


def _check_preserves_subgroup(automorphism: Automorphism, subgroup: Subgroup) -> bool:
    if subgroup is not None:
        return all((automorphism(h) in subgroup for h in subgroup.generators))
    return True


def _check_order_condition(group: Group, in_out_dict: dict) -> bool:
    for a, b in product(group.generators, group.generators):
        order_condition = (a * b).order == (in_out_dict[a] * in_out_dict[b]).order
        if not order_condition:
            return False
    return True


def _check_bijection_condition(group: Group, generator_in_out_dict: dict) -> bool:
    image = group.subgroup_generated_by(list(generator_in_out_dict.values()))
    if image.order != group.order:
        return False
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


def _fetch_potential_automorphisms(group: Group) -> list:
    possible_generator_images = {
        g: [h for h in group if g.order == h.order] for g in group.generators
    }
    generator_image_choices = product(
        *[possible_generator_images[g] for g in possible_generator_images.keys()]
    )
    return [
        {g: image_list[i] for i, g in enumerate(possible_generator_images.keys())}
        for image_list in generator_image_choices
    ]


def _filter_potential_automorphisms(
    automorphisms: list, potential_automorphisms: list, group: Group
) -> list:
    in_out_dicts_to_remove = [
        {g: inn(g) for g in group.generators} for inn in automorphisms
    ]
    potential_automorphisms = [
        in_out_dict
        for in_out_dict in potential_automorphisms
        if in_out_dict not in in_out_dicts_to_remove
    ]
    potential_automorphisms = [
        in_out_dict
        for in_out_dict in potential_automorphisms
        if _check_order_condition(group, in_out_dict)
    ]
    potential_automorphisms = [
        in_out_dict
        for in_out_dict in potential_automorphisms
        if _check_bijection_condition(group, in_out_dict)
    ]

    return potential_automorphisms


def Inn(group: Group, relative_subgroup: Subgroup = None) -> Group:
    """Creates the group Inn(G) of inner automorphisms of the finite group G."""
    inner_automorphisms = list(
        {
            Automorphism(group, _inner_automorphism_factory(g))
            for g in [group.identity] + [x for x in group if x not in group.center]
        }
    )
    if relative_subgroup is not None:
        inner_automorphisms = [
            inner_auto
            for inner_auto in inner_automorphisms
            if _check_preserves_subgroup(inner_auto, relative_subgroup)
        ]
    inner_automorphism_group = Group(inner_automorphisms)
    inner_automorphism_group.identity = _aut_get_identity(group)
    return inner_automorphism_group


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
    automorphisms = Inn(group, relative_subgroup).elements

    potential_automorphisms = _fetch_potential_automorphisms(group)
    potential_automorphisms = _filter_potential_automorphisms(
        automorphisms, potential_automorphisms, group
    )
    updated_potential_automorphisms = copy(potential_automorphisms)

    for in_out_dict in potential_automorphisms:
        if in_out_dict in updated_potential_automorphisms:
            automorphism = Automorphism.from_action_on_generators(group, in_out_dict)

            homomorphism_check = automorphism.validate_homomorphism()
            preserves_conjugacy_classes = _check_class_preserving(
                automorphism, class_preserving
            )
            fixes_subgroup = _check_preserves_subgroup(automorphism, relative_subgroup)

            if all(
                [
                    homomorphism_check,
                    preserves_conjugacy_classes,
                    fixes_subgroup,
                ]
            ):
                automorphisms = list(
                    {A * automorphism for A in automorphisms}.union(set(automorphisms))
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


def Out(group: Group, relative_subgroup: Subgroup) -> Group:
    """
    Creates the group Out(G) of outer automorphisms of the finite group G.

    Out(G) is defined as the quoteint group Aut(G)/Inn(G).
    """
    automorphism_group = Aut(group, relative_subgroup)
    inner_automorphism_group = Inn(group, relative_subgroup)
    return automorphism_group / Subgroup(
        inner_automorphism_group.elements, automorphism_group
    )
