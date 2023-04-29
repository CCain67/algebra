import pytest

from group_theory.common_groups.residue_class_groups import *
from group_theory.common_groups.permutation_groups import *
from group_theory.common_groups.matrix_groups import *
from group_theory.common_groups.misc_groups import *


def test_symmetric_quotient_alternating():
    S = symmetric_group(4)
    A = S.commutator_subgroup
    assert (S/A).order == 2