from itertools import product
from math import lcm

from group_theory.base.group_elements import (
    GroupElement,
)
from group_theory.base.groups import Group

class CartesianProductElement(GroupElement):
    def __init__(self, elements: tuple[GroupElement]) -> None:
        self.elements = elements
        self.num_elements = len(self.elements)

    def __repr__(self):
        return str(tuple(x for x in self.elements))
    
    def __hash__(self):
        return hash(self.elements)
    
    def __eq__(self,other):
        return self.elements==other.elements

    def __mul__(self, other):
        if len(self.elements)!=len(other.elements):
            raise ValueError('the length of the cartesian products do not match')
        product = tuple(x*y for x,y in zip(self.elements, other.elements))
        return CartesianProductElement(product)
    
    def __invert__(self):
        inv = tuple(~x for x in self.elements)
        return CartesianProductElement(inv)
    
    def __getitem__(self, key): 
        return self.elements[key]
    
    def __iter__(self): 
        return iter(self.elements)

    def get_order(self):
        return lcm(*[x.order for x in self.elements])
    
    def is_identity(self):
        return [x.is_identity() for x in self.elements] == [True]*self.num_elements
    
class CartesianProduct(Group):
    def __init__(self, *args: Group):
        self.list_of_groups = args
        self.elements = self.get_elements()
        self.identity = self.get_identity()
        self.order = len(self.elements)
        
    def get_elements(self):
        group_elements = [G.elements for G in self.list_of_groups ]
        product_elements = list(CartesianProductElement(g) for g in product(*group_elements))
        return product_elements