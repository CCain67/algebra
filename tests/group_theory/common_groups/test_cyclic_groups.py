"""Unit tests for cyclic groups"""

from noetherpy.group_theory.common_groups.cyclic_groups import (
    cyclic_group,
)


def test_cyclic_group_order():
    assert cyclic_group(31).order == 31
