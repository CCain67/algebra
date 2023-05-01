import pytest

from group_theory.group_elements import Permutation


class Base:
    permutation = None

    def setup_method(self):
        self.P = Permutation(self.permutation)

    def test_order(self, benchmark):
        benchmark(self.P.get_order)
    
    def test_cycle_decomp(self, benchmark):
        benchmark(self.P.get_cycle_decomposition)

    def test_cycle_notation(self, benchmark):
        benchmark(self.P.get_cycle_notation)


@pytest.mark.benchmark(group="Permutation: product of transpositions on 10 letters")
class Test_Permutation_10_letter_transpositions(Base):
    permutation = {1:2, 2:1, 3:4, 4:3, 5:6, 6:5, 7:8, 8:7, 9:10, 10:9}

@pytest.mark.benchmark(group="Permutation: cycle on 10 letters")
class Test_Permutation_10_letter_cycle(Base):
    permutation = {1:2, 2:3, 3:4, 4:5, 5:8, 6:7, 7:9, 8:6, 9:10, 10:1}

@pytest.mark.benchmark(group="Permutation: cycle times transpositions on 10 letters")
class Test_Permutation_10_letter_cycle_times_transpositions(Base):
    permutation = {1:2, 2:1, 3:4, 4:3, 5:6, 6:5, 7:8, 8:9, 9:10, 10:7}