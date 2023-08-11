"""
This module defines functions which compute centralizers and normalizers of subgroups. 
"""
from typing import (
    Iterable,
)
from ..group_elements import GroupElement
from ..groups import (
    Group,
    Subgroup,
)


def centralizer(subgroup: Subgroup, group: Group) -> Subgroup:
    """Computes the centralizer of a subgroup.

    Args:
        H (Subgroup): A subgroup of of the group (group)

    Raises:
        ValueError: if the parent group of the subgroup is not (group), then
        this error is raised.

    Returns:
        Subgroup: the centralizer subgroup of the subgroup passed.
    """
    if group != subgroup.parent_group:
        raise ValueError("subgroup provided is not a subgroup of this group")
    centralizer_elements = [group.identity]
    for z in (x for x in group if x != group.identity):
        if all((z * h == h * z for h in subgroup.generators)):
            centralizer_elements += [z]
    return Subgroup(centralizer_elements, group)


def normalizer(subgroup: Subgroup, group: Group) -> Subgroup:
    """Computes the normalizer of a subgroup.

    Args:
        subgroup (Subgroup): the subgroup whose normalizer we want to compute.

    Raises:
        ValueError: if the parent group of the subgroup is not (group), then
        this error is raised.

    Returns:
        Subgroup: the normalizer of the subgroup provided.
    """
    if group != subgroup.parent_group:
        raise ValueError("subgroup provided is not a subgroup of this group")
    normalizer_elements = [group.identity]
    for z in (x for x in group if x != group.identity):
        if all((z * h * (~z) in subgroup for h in subgroup.generators)):
            normalizer_elements += [z]
    return Subgroup(normalizer_elements, group)


def conjugate_subgroup(subgroup: Subgroup, g: GroupElement) -> Subgroup:
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
    if g not in subgroup.parent_group:
        raise ValueError("group element must be a member of the parent group")
    return subgroup.parent_group.subgroup_generated_by(
        [g * h * (~g) for h in subgroup.generators]
    )


def normal_closure(elements: Iterable, group: Group) -> Subgroup:
    """Computes the normal closure of a subset of a group.

    Args:
        elements (Iterable): the set of elements to compute the normal closure of.
        group (Group): the parent group of the set of elements

    Returns:
        Subgroup: the normal closure of the set of elements provided.
    """
    return group.subgroup_generated_by(
        [g * s * (~g) for s in elements for g in group], group
    )


def normal_core(subgroup: Subgroup) -> Subgroup:
    """Computes the normal core of the subgroup provided. The normal core of a subgroup
    is defined as the intersection of all conjugates of the subgroup. It suffices to
    compute this only on generators of the group.

    Args:
        subgroup (Subgroup): the subgroup to compute the normal core of

    Returns:
        Subgroup: the normal core of the subgroup
    """
    core = subgroup
    for g in subgroup.parent_group.generators:
        core = core & conjugate_subgroup(subgroup, g)
    return core
