import pytest

from group_theory.common_groups.permutation_groups import SymmetricGroup


class Base:
    order = None

    def setup_method(self):
        self.Sym = SymmetricGroup.as_permutation_group(self.order)

    def test_load(self, benchmark):
        benchmark(self.Sym.as_permutation_group,self.order)

    def test_center(self, benchmark):
        benchmark(self.Sym.get_center)


@pytest.mark.benchmark(group="SymmetricGroup on 4 letters")
class Test_SymmetricGroup_4(Base):
    order = 4

@pytest.mark.benchmark(group="SymmetricGroup on 5 letters")
class Test_SymmetricGroup_5(Base):
    order = 5