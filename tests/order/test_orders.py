import pytest

from group_theory.common_groups.residue_class_groups import *
from group_theory.common_groups.permutation_groups import *
from group_theory.common_groups.matrix_groups import *
from group_theory.common_groups.misc_groups import *

def test_cyclic_group_order():
    assert CyclicGroup.as_residue_class_group(29).order == 29