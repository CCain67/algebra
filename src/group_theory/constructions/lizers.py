"""
This module defines functions which compute centralizers and normalizers of subgroups. 
"""

from group_theory.groups import (
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
    number_of_generators = len(subgroup.generators)
    for z in [x for x in group if x != group.identity]:
        counter = 0
        for h in subgroup.generators:
            if z * h == h * z:
                counter += 1
                continue
            else:
                break
        if counter == number_of_generators:
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
    for z in [x for x in group if x != group.identity]:
        if subgroup.left_coset(z) == subgroup.right_coset(z):
            normalizer_elements += [z]
    return Subgroup(normalizer_elements, group)
