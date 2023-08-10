"""Unit tests for misc groups (e.g. polyhedral groups, Klein 4 group, etc.)"""

from noetherpy.group_theory.common_groups.misc_groups import (
    dihedral_group,
)


def test_dihedral_group():
    assert dihedral_group(7).order == 14
