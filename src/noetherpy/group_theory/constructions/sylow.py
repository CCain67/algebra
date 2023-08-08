"""This module contains functions for performing basic Sylow theory."""

from galois import factors

from ..groups import (
    Group,
    Subgroup,
)
from ..constructions.subgroup_constructions import (
    centralizer,
    conjugate_subgroup,
    normalizer,
)


def sylow_p_subgroup(group: Group, p: int) -> Subgroup:
    """Computes a random Sylow p-subgroup of the given group.

    Args:
        group (Group): The parent group.
        p (int): prime dividing the order of the group given.

    Raises:
        ValueError: This error is raised if the prime given does not
        divide the order of the group.

    Returns:
        Subgroup: A Sylow p-subgroup of the group given.
    """
    prime_factors = factors(group.order)
    prime_decomp = dict(list(zip(prime_factors[0], prime_factors[1])))
    if p not in prime_decomp:
        raise ValueError("the order of the group must be divisible by the prime given")
    required_subgroup_order = p ** prime_decomp[p]

    exp, exp_max = 1, prime_decomp[p]

    sylow_p_gens = [group.element_by_order(p)]
    sylow_subgroup = group.subgroup_generated_by(sylow_p_gens)

    while (sylow_subgroup.order != required_subgroup_order) and (exp <= exp_max):
        found_element = False
        for g in normalizer(sylow_subgroup, group):
            if (g not in sylow_subgroup) and (g.order == p**exp):
                sylow_p_gens.append(g)
                found_element = True
                break
        if not found_element:
            exp += 1
            sylow_p_gens = [group.element_by_order(p**exp)]
        sylow_subgroup = group.subgroup_generated_by(sylow_p_gens)

    assert sylow_subgroup.order == required_subgroup_order

    return sylow_subgroup


def number_of_sylow_p_subgroups(group: Group, p: int) -> int:
    """Computes the number of Sylow p-subgroups of the given group

    Args:
        group (Group): the group to compute the number of Sylow p-subgroups for
        p (int): the prime dividing the order of the group given.

    Returns:
        int: the number of Sylow p-subgroups of the given group.
    """
    return int(group.order / normalizer(sylow_p_subgroup(group, p), group).order)


def fetch_all_sylow_p_subgroups(group: Group, p: int) -> list:
    """Generates the list of all Sylow p-subgroups of the group given.

    Args:
        group (Group): The group whose Sylow p-subgroups we want to compute
        p (int): prime number dividing the order of the group.

    Raises:
        ValueError: This error is raised if the prime given does not
        divide the order of the group.

    Returns:
        list: list of all Sylow p-subgroups of the given group.
    """
    prime_factors = factors(group.order)
    prime_decomp = dict(list(zip(prime_factors[0], prime_factors[1])))
    if p not in prime_decomp:
        raise ValueError("the order of the group must be divisible by the prime given")

    sylow_p = sylow_p_subgroup(group, p)
    sylow_centralizer = centralizer(sylow_p, group)
    sylow_p_subgroups = set()
    for g in (x for x in group if x not in sylow_centralizer):
        sylow_p_subgroups.update({conjugate_subgroup(sylow_p, g)})

    return list(sylow_p_subgroups)
