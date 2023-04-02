from typing import Callable, Tuple

from group_theory.base.group_elements import (
    CartesianProductElement,
    GroupElement,
)
from group_theory.base.groups import (
    Group,
    Subgroup,
)

class Homomorphism:
    def __init__(self, domain: Group, map: Callable[[Tuple[GroupElement,...]],GroupElement], codomain: Group) -> None:
        self.domain = domain
        self.map = map
        self.codomain = codomain

        # properties
        self._graph = None
        self._image = None
        self._kernel = None
        self._is_iso = None

    def validate_homomorphism(self) -> bool:
        return self.graph.order==self.domain.order
    
    @classmethod
    def from_dict(cls, domain: Group, in_out_dict: dict, codomain: Group):
        def f(g: GroupElement) -> GroupElement:
            return in_out_dict[g]
        
        return cls(domain, f, codomain)
    
    def get_graph(self) -> Subgroup:
        self.domain.get_random_generators()
        augmented_generators = self.domain.generators+[a*b for a in self.domain.generators for b in self.domain.generators]
        G_times_H = self.domain*self.codomain
        element_image_pairs = [CartesianProductElement((x,self.map(x))) for x in augmented_generators]
        X = G_times_H.subgroup_generated_by(element_image_pairs)
        return X

    @property
    def graph(self):
        if self._graph is None:
            self._graph = self.get_graph()
            return self._graph
        else:
            return self._graph
    
    def get_image(self) -> Group:
        return NotImplemented

    @property
    def image(self):
        if self._image is None:
            self._image = self.get_image()
            return self._image
        else:
            return self._image
        
    def get_kernel(self) -> Subgroup:
        X = self.graph
        G_times_H = self.domain*self.codomain
        G_1 = Subgroup([p for p in G_times_H if p[1].is_identity()], G_times_H)
        return Subgroup([p[0] for p in (X & G_1)], self.domain)

    @property
    def kernel(self):
        if self._kernel is None:
            self._kernel = self.get_kernel()
            return self._kernel
        else:
            return self._kernel
        
    def check_iso(self) -> bool:
        '''
        this utilizes the fact that for two finite sets, if an injective map exists between them, and 
        the two sets are of equal order, then the map is a bijection. 
        '''
        if len(self.kernel)==1 and self.domain.order==self.codomain.order and self.validate_homomorphism():
            return True
        else:
            return False

    @property
    def is_iso(self):
        if self._is_iso is None:
            self._is_iso = self.check_iso()
            return self._is_iso
        else:
            return self._is_iso