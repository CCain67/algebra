from __future__ import annotations

from random import sample
from itertools import product

from group_theory.base.group_elements import (
    CartesianProductElement,
    GroupElement,
)

class Group:
    def __init__(self, elements: list[GroupElement]):
        self.elements = elements
        self.order = len(self.elements)
        self.identity = self.get_identity()
        self.canonical_generators = None

        # properties
        self._generators = None
        self._center = None
        self._commutator_subgroup = None
        self._conjugacy_classes = None

    def get_identity(self) -> GroupElement:
        for g in self.elements:
            if g.is_identity():
                self.identity = g
                return g

    def __repr__(self):
        s = ""
        for g in self:
            s+=repr(g)+'\n'
        return s

    def __mul__(self,other):
        '''
        this is the standard cartesian product of groups: for group G and H, G*H contains all
        pairs of elements of the form (g,h)
        '''
        product_elements = list(CartesianProductElement(g) for g in product(self.elements,other.elements))
        return Group(product_elements)

    def __truediv__(self,other):
        if isinstance(other, Subgroup):
            if other.parent_group==self:
                cosets = {Coset(g, other) for g in self}
                return Group(list(cosets))
            else:
                raise ValueError('the subgroup provided is not a subgroup of the provided parent group')
        else:
            raise TypeError('you must pass a valid subgroup to quotient by')

    def __eq__(self, other):
        return set(self.elements)==set(other.elements)
    
    def __ne__(self, other):
        return set(self.elements)!=set(other.elements)
    
    def __lt__(self, other):
        return set(self.elements).issubset(set(other.elements))
    
    def __le__(self, other):
        return set(self.elements).issubset(set(other.elements))
    
    def __gt__(self, other):
        return set(self.elements).issuperset(set(other.elements))
    
    def __ge__(self, other):
        return set(self.elements).issuperset(set(other.elements))
    
    def __hash__(self):
        return hash(frozenset(self.elements))
    
    def __getitem__(self, key): 
        return self.elements[key]
    
    def __iter__(self): 
        return iter(self.elements)
    
    def __len__(self):
        return len(self.elements)
    
    def is_trivial(self):
        if len(self.elements)==1 and self.elements[0].is_identity():
            return True
        else:
            return False
        
    def is_abelian(self):
        return self.commutator_subgroup.is_trivial()
    
    def is_perfect(self):
        return self.commutator_subgroup == self
    
    def is_solvable(self):
        return self.derived_series()[-1].is_trivial()
    
    def is_nilpotent(self):
        return self.lower_central_series()[-1].is_trivial()
    
    def get_random_generators(self) -> list[GroupElement]:
        if len(self.elements)==1:
            return [self[0]]
        generators = []
        elements = [g for g in self if not g.is_identity()]
        generators+=[sample(elements,1)[0]]
        generated_elements = self.subgroup_generated_by(generators).elements
        while len(generated_elements)!=self.order:
                elements = [g for g in elements if g not in generated_elements]
                generators+=[sample(elements,1)[0]]
                generated_elements = self.subgroup_generated_by(generators).elements
        return generators
    
    def get_generators(self):
        if self.canonical_generators:
            return self.canonical_generators
        else:
            return self.get_random_generators()

    @property
    def generators(self):
        if self._generators is None:
            self._generators = self.get_generators()
            return self._generators
        else:
            return self._generators
    
    def get_center(self):
        '''
        the centralizer of the entire group is the center
        '''
        G = Subgroup(self.elements,self)
        return self.centralizer(G)
    
    @property
    def center(self):
        if self._center is None:
            self._center = self.get_center()
            return self._center
        else:
            return self._center

    def get_conjugacy_classes(self):
        '''
        each element of the center forms a conjugacy class consisting of just itself, so we compute the center
        first and the other classes second
        '''
        conjugacy_classes = set()
        ZG = self.center
        conjugacy_classes.update([frozenset({z}) for z in ZG])
        for h in set(self.elements)-set(ZG.elements):
            conjugacy_classes.add(frozenset({g*h*(~g) for g in self}))
        return conjugacy_classes
    
    @property
    def conjugacy_classes(self):
        if self._conjugacy_classes is None:
            self._conjugacy_classes = self.get_conjugacy_classes()
            return self._conjugacy_classes
        else:
            return self._conjugacy_classes
    
    def get_commutator_subgroup(self):
        if self.generators:
            '''
            if G is generated by a set S, then the commutator subgroup is the normal closure of the set of commutators of elements of S
            '''
            return self.subgroup_generated_by([g*x*y*(~x)*(~y)*(~g) for x in self.generators for y in self.generators for g in self.elements])
        else:
            return self.subgroup_generated_by([x*y*(~x)*(~y) for x in self.elements for y in self.elements])

    @property
    def commutator_subgroup(self):
        if self._commutator_subgroup is None:
            self._commutator_subgroup = self.get_commutator_subgroup()
            return self._commutator_subgroup
        else:
            return self._commutator_subgroup
    
    def centralizer(self, H: Subgroup):
        if self!=H.parent_group:
            raise ValueError('subgroup provided is not a subgroup of this group')
        Z = [self.identity]
        number_of_generators = len(H.generators)
        for z in [x for x in self if x!=self.identity]:
            counter = 0 
            for h in H.generators:
                if z*h==h*z:
                    counter += 1
                    continue
                else:
                    break
            if counter == number_of_generators:
                Z += [z]
        return Subgroup(Z,self)
    
    def normalizer(self, H: Subgroup):
        if self!=H.parent_group:
            raise ValueError('subgroup provided is not a subgroup of this group')
        Z = [self.identity]
        for z in [x for x in self if x!=self.identity]:
            if H.left_coset(z)==H.right_coset(z):
                Z+=[z]
        return Subgroup(Z,self)

    def subgroup_generated_by(self, generators: list[GroupElement]) -> Subgroup:
        if type(generators)!=list:
            raise TypeError('the generators must be a list of group elements')
        '''
        this algorithm is the "black-box" algo found here:
        https://groupprops.subwiki.org/w/index.php?title=Black-box_group_algorithm_for_finding_the_subgroup_generated_by_a_subset
        '''
        e = self.identity
        H = {e}
        F = {e}
        
        while F:
            K = {f*s for f in F for s in generators}
            K = K-H
            F = K
            H = F.union(H)
            
        return Subgroup(list(H), self)
    
    def derived_series(self):
        last_subgroup = Subgroup(self.elements, self)
        derived_series_list = [last_subgroup]
        next_subgroup = last_subgroup.commutator_subgroup
        while next_subgroup != last_subgroup:
            derived_series_list.append(next_subgroup)
            last_subgroup = next_subgroup
            next_subgroup = last_subgroup.commutator_subgroup
        return derived_series_list

    def lower_central_series(self):
        last_subgroup = Subgroup(self.elements, self)
        lower_central_series_list = [last_subgroup]
        next_subgroup = last_subgroup.commutator_subgroup
        while next_subgroup != last_subgroup:
            lower_central_series_list.append(next_subgroup)
            last_subgroup = next_subgroup
            next_subgroup = self.subgroup_generated_by([x*y*(~x)*(~y) for x in last_subgroup for y in self])
        return lower_central_series_list

    def upper_central_series(self):
        last_center = Subgroup([self.identity], self)
        upper_central_series_list = [last_center]
        next_center = self.center()
        while last_center != next_center:
            upper_central_series_list.append(next_center)
            last_center = next_center
            next_center = Subgroup(
                [x for x in self if {x*y*(~x)*(~y) for y in self}.issubset(last_center)], self
            )
        return upper_central_series_list

class Subgroup(Group):
    def __init__(self, elements: list[GroupElement], parent_group: Group):
        super().__init__(elements)
        self.parent_group = parent_group
        self.canonical_generators = None
        self.index = int(self.parent_group.order/self.order)

        # properties
        self._is_normal = None

    def validate_inclusion(self):
        return set(self.elements).issubset(self.parent_group.elements)
    
    # IMPORTANT - until this is called and returns true, 
    # the Subgroup instance may not actually be a subgroup
    def validate_subgroup(self):
        for g in self:
            for h in self:
                if g*(~h) not in self:
                    return False
        return True
    
    def __repr__(self):
        return str(self.elements)
    
    def __eq__(self, other):
        return set(self.elements)==set(other.elements)
    
    def __hash__(self):
        return hash((self.elements,self.parent_group))
    
    def __and__(self, other):
        if self.parent_group!= other.parent_group:
            raise ValueError('the subgroups must both be a subset of the same group')
        return Subgroup([x for x in self.elements if x in other.elements], self.parent_group)
    
    def __matmul__(self, other):
        if self.parent_group!= other.parent_group:
            raise ValueError('the subgroups must both be a subset of the same group')
        HK = {h*k for h in self for k in other}
        KH = {k*h for h in self for k in other}
        if HK==KH:
            return Subgroup(list(HK), self.parent_group)
        else: 
            return HK
    
    def left_coset(self, g: GroupElement):
        if g not in self.parent_group:
            raise ValueError('group element must be a member of the parent group')
        return {g*x for x in self.elements}
    
    def right_coset(self, g: GroupElement):
        if g not in self.parent_group:
            raise ValueError('group element must be a member of the parent group')
        return {x*g for x in self.elements}
    
    def conjugate_subgroup(self, g: GroupElement):
        if g not in self.parent_group:
            raise ValueError('group element must be a member of the parent group')
        return Subgroup([g*x*(~g) for x in self], self.parent_group)
    
    def check_normality(self):
        if self.index == 2: # index 2 subgroups are always normal
            return True
        for g in self.parent_group.generators:
            if self.left_coset(g)!=self.right_coset(g):
                return False
        return True

    @property
    def is_normal(self):
        if self._is_normal is None:
            self._is_normal = self.check_normality()
            return self._is_normal
        else:
            return self._is_normal

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