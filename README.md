# NoetherPy
A small python library for computational abstract algebra designed to make it easy to learn, test conjectures and ideas in a python notebook. 

## Features
Most of basic (finite) group theory has been implemented: 
- Simple and intuitive ways to define groups, subgroups, homomorphisms, etc. 
- Symmetric/alternating groups, matrix groups over finite fields, cyclic groups and groups of units of (implemented) rings, dihedral groups and the Klein 4 group are all implemented.
  - Many of the above groups can be instantiated as matrix groups, or permutation groups.
- Most well-known operations on groups/subgroups are supported:
  - Cartesian products, quotient groups, and semidirect products
  - centers, centralizers, normalizers, normal cores, etc.
  - basic Sylow theory 
  - automorphism groups (along with inner and outer automorphism groups)
  - subgroup series (derived, lower and upper central)
- Group actions can be defined on anything of Iterable type
- Utilizes the [galois](https://github.com/mhostetter/galois) package for finite field / number theory support 

## Future plans:
- better documentation + wiki
- representation and character theory of finite groups
- group (co)homology

#### You know GAP/Magma/Sage/etc. exists right??
Of course! Those tools are certainly better/faster, this is just for fun.

