from group_theory.base.group_elements import GroupElement
from group_theory.base.groups import (
    Group,
    Subgroup,
)


class Coset(GroupElement):
    '''
    this class is used to represent elements gH of quotient groups G/H 
    for elements g and normal subgroups H
    '''
    def __init__(self, g: GroupElement, subgroup: Subgroup) -> None:
        self.g, self.subgroup = self.validate_coset(g, subgroup)
        self.elements = [self.g*x for x in self.subgroup]

    def __repr__(self):
        return str(self.g)+'H'
    
    def __eq__(self, other):
        return set(self.elements)==set(other.elements)
    
    def __ne__(self, other):
        return set(self.elements)!=set(other.elements)
    
    def __hash__(self):
        return hash(frozenset(self.elements))
    
    def __getitem__(self, key): 
        return self.elements[key]
    
    def __iter__(self): 
        return iter(self.elements)
    
    def __len__(self):
        return len(self.elements)
    
    def __mul__(self,other):
        if self.subgroup!=other.subgroup:
            raise ValueError('the subgroups of each coset must match')
        return Coset(self.g*other.g, self.subgroup)
    
    def __invert__(self):
        return Coset(~self.g, self.subgroup)

    @staticmethod
    def validate_coset(g: GroupElement, subgroup: Subgroup):
        assert g in subgroup.parent_group, "the group element must lie in the parent group"
        return g, subgroup

    def get_order(self):
        return NotImplemented
    
    def is_identity(self):
        return self.g in self.subgroup
    
class QuotientGroup(Group):
    def __init__(self, parent_group: Group, normal_subgroup: Subgroup):
        self.parent_group = parent_group
        self.normal_subgroup = normal_subgroup
        self.elements = self.get_elements()
        self.identity = self.get_identity()
        self.order = len(self.elements)
        
    def get_elements(self):
        cosets = {Coset(g, self.normal_subgroup) for g in self.parent_group}
        return list(cosets)