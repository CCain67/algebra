# NoetherPy
A python library for computational abstract algebra designed to make it easy to learn, test conjectures and ideas in a python notebook. 

<div align=center>
  <a href="https://github.com/CCain67/noetherpy/actions/workflows/test.yml"><img src="https://github.com/CCain67/noetherpy/actions/workflows/test.yml/badge.svg"></a>
  <a href="https://github.com/CCain67/noetherpy/actions/workflows/pylint.yml"><img src="https://github.com/CCain67/noetherpy/actions/workflows/pylint.yml/badge.svg"></a>
</div>

## Features
Most of basic (finite) group theory has been implemented: 
- Simple and intuitive ways to define groups, subgroups, homomorphisms, etc. 
- Symmetric/alternating groups, matrix groups over finite fields, cyclic groups and groups of units of (implemented) rings, finite abelian groups, dihedral groups, quasidihedral groups, and the Klein 4 group are all implemented.
  - Many of the above groups can be instantiated as matrix groups, or permutation groups.
- Most well-known operations on groups/subgroups are supported:
  - Cartesian products, quotient groups, and semidirect products
  - centers, centralizers, normalizers, normal cores, etc.
  - basic Sylow theory 
  - automorphism groups (along with inner and outer automorphism groups)
  - subgroup series (derived, lower and upper central)
- Group actions can be defined on anything of Iterable type
- Utilizes the [galois](https://github.com/mhostetter/galois) package for finite field / number theory support 

## Getting Started

First, clone this repo:
```
git clone https://github.com/CCain67/noetherpy.git
```
then, for now it is recommended to work in a python 3.11 virtual environment:
```
conda create -n venv_name python=3.11
```
and install the package by running the following command in the project's root directory:
```
pip install .
```
Now, inside of a Jupyter Notebook, the package can be imported:
```python
import noetherpy as ntr
```

### Showcase

For a getting a general idea of what one can do with this library, check the showcase here:
[Showcase](showcase/)

## Future plans:
- better documentation + wiki
- representation and character theory of finite groups
- group (co)homology

#### You know GAP/Magma/Sage/etc. exists right??
Of course! Those tools are certainly better/faster, this is just for fun.

