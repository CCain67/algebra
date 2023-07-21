"""
This module defines functions which compute various subgroup series 
of a group.
"""

from group_theory.groups import (
    Group,
    Subgroup,
)


def derived_series(group: Group) -> list[Subgroup]:
    """Computes the derived series of the group.

    Returns:
        list: a list of Subgroups which make up the derived series of the group.
    """
    last_subgroup = Subgroup(group.elements, group)
    derived_series_list = [last_subgroup]
    next_subgroup = last_subgroup.commutator_subgroup
    while next_subgroup != last_subgroup:
        derived_series_list.append(next_subgroup)
        last_subgroup = next_subgroup
        next_subgroup = last_subgroup.commutator_subgroup
    return derived_series_list


def lower_central_series(group: Group) -> list[Subgroup]:
    """Computes the lower central series of the group.

    Returns:
        list: a list of Subgroups which make up the lower central series of the group.
    """
    last_subgroup = Subgroup(group.elements, group)
    lower_central_series_list = [last_subgroup]
    next_subgroup = last_subgroup.commutator_subgroup
    while next_subgroup != last_subgroup:
        lower_central_series_list.append(next_subgroup)
        last_subgroup = next_subgroup
        next_subgroup = group.subgroup_generated_by(
            [x * y * (~x) * (~y) for x in last_subgroup for y in group]
        )
    return lower_central_series_list


def upper_central_series(group: Group) -> list[Subgroup]:
    """Computes the upper central series of the group.

    Returns:
        list: a list of Subgroups which make up the upper central series of the group.
    """
    last_center = Subgroup([group.identity], group)
    upper_central_series_list = [last_center]
    next_center = group.center
    while last_center != next_center:
        upper_central_series_list.append(next_center)
        last_center = next_center
        next_center = Subgroup(
            [
                x
                for x in group
                if {x * y * (~x) * (~y) for y in group}.issubset(last_center)
            ],
            group,
        )
    return upper_central_series_list


def is_solvable(group) -> bool:
    """Determines whether or not the group is solvable, i.e., if the derived
    series of G(group) ends at the trivial group.

    Returns:
        bool: True if the group is solvable, False otherwise.
    """
    if group.order % 2 == 1:  # the is the Feit-Thompson Theorem
        return True
    return derived_series(group)[-1].is_trivial()


def is_nilpotent(group) -> bool:
    """Determines whether or not the group is nilpotent, i.e., if the lower
    central series of G(group) ends at the trivial group.

    Returns:
        bool: True if the group is trivial, False otherwise.
    """
    return lower_central_series(group)[-1].is_trivial()
