"""Unit tests for Sylow theory functions"""

import pytest

from noetherpy.group_theory.constructions.sylow import (
    random_sylow_p_subgroup,
    number_of_sylow_p_subgroups,
    fetch_all_sylow_p_subgroups,
)

from noetherpy.group_theory.common_groups.permutation_groups import (
    alternating_group,
    symmetric_group,
)
from noetherpy.group_theory.common_groups.misc_groups import (
    dihedral_group,
    klein_four_group,
)


dih_7 = dihedral_group(7)
dih_8 = dihedral_group(8)

sym_5 = symmetric_group(5)
alt_4 = alternating_group(4)
klein_4 = klein_four_group(representation="permutation")


def test_random_sylow_p_subgroup() -> None:
    # the alternating group A_4 has a unique Sylow 2-subgroup,
    # so the randomness here is not an issue
    assert random_sylow_p_subgroup(alt_4, 2) == klein_4


def test_number_of_sylow_p_subgroups_computation_gives_correct_value() -> None:
    assert number_of_sylow_p_subgroups(dih_7, 2) == 7
    assert number_of_sylow_p_subgroups(sym_5, 5) == 6
    assert number_of_sylow_p_subgroups(alt_4, 3) == 4


def test_fetch_all_sylow_p_subgroups_finds_correct_number_of_sylow_subgroups() -> None:
    assert len(fetch_all_sylow_p_subgroups(dih_7, 2)) == 7
    assert len(fetch_all_sylow_p_subgroups(sym_5, 5)) == 6
    assert len(fetch_all_sylow_p_subgroups(alt_4, 3)) == 4
