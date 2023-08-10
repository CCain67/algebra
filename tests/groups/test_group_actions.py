"""Unit tests for group actions"""

from noetherpy.group_theory.group_elements import Permutation
from noetherpy.group_theory.group_actions import GroupAction
from noetherpy.group_theory.common_groups.permutation_groups import symmetric_group


def test_symmetric_group_action() -> None:
    G = symmetric_group(4)
    X = {0, 1, 2, 3}

    def action(g: Permutation, x: int) -> int:
        return g[x]

    sym_4_action = GroupAction(G, X, action)

    assert sym_4_action.validate_action()
    assert sym_4_action.is_transitive
    assert sym_4_action.is_faithful
    assert not sym_4_action.is_free
