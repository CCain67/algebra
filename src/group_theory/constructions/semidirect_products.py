from itertools import product

from group_theory.groups import (
    Group,
)
from group_theory.group_elements import GroupElement
from group_theory.homomorphisms import Homomorphism

class SemidirectProductElement(GroupElement):
    '''
    Here, we are viewing the 'twisting' homomorphism twist:G -> Aut(H) as a map
    twist:G x H -> H which restricts to an automorphism of H for each g in G.
    '''
    def __init__(self, elements: tuple[GroupElement], twist: Homomorphism) -> None:
        self.elements = elements
        self.twist = twist
        self.num_elements = 2

        # properties
        self._order = None

    def __repr__(self):
        return str(self.elements)
    
    def __hash__(self):
        return hash((self.elements,self.twist))
    
    def __eq__(self,other):
        return self.elements==other.elements
    
    def __ne__(self,other):
        return self.elements!=other.elements

    def __mul__(self, other):
        if self.twist != other.twist:
            raise ValueError('the twisting homomorphisms do not match')
        product = (self.elements[0]*self.twist(self.elements[1],other.elements[0]), self.elements[1]*other.elements[1])
        return SemidirectProductElement(product, self.twist)
    
    def __invert__(self):
        return SemidirectProductElement((self.twist(~self.elements[1],~self.elements[0]), ~self.elements[1]), self.twist)
    
    def __getitem__(self, key): 
        return self.elements[key]
    
    def __iter__(self): 
        return iter(self.elements)

    def get_order(self) -> int:
        A = self
        i=1
        while not A.is_identity():
            A *= self
            i += 1
        return i
    
    @property
    def order(self) -> int:
        if self._order is None:
            self._order = self.get_order()
            return self._order
        else:
            return self._order
    
    def is_identity(self):
        return self.elements[0].is_identity() and self.elements[1].is_identity()
    

class SemidirectProduct(Group):
    '''
    Terminology: for a group G with action f:G x H -> H, we will refer to G as the "base" group and H as the "fiber" group, in an analogy to fiber bundles.
    '''
    def __init__(self, fiber: Group, twist: Homomorphism, base: Group):
        self.base = base
        self.fiber = fiber
        self.twist = twist
        
        elements = list(SemidirectProductElement(p,self.twist) for p in product(self.fiber.elements,self.base.elements))
        super().__init__(elements)

