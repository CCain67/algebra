import pytest 

from group_theory.common_groups.misc_groups import DihedralGroup

def test_center_dihedral():
    assert DihedralGroup.as_permutation_group(7).center().order == 1
    assert DihedralGroup.as_permutation_group(8).center().order == 2