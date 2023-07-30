# The Exotic Outer Automorphism of $S_6$

```python
import noetherpy as ntr
```

It is well known that the 6th symmetric group $S_6$ is the only symmetric group which has a non-trivial outer automorphism group. In particular, $\mathrm{Out}(S_6) \cong \mathbb{Z}/2\mathbb{Z}$. This is easily verified with `noetherpy` by calling:

```python
G = ntr.symmetric_group(6)
OutG = ntr.Out(G)
```

From here, calling `OutG` shows two elements, as expected:

```python
Coset(
(0 1 2 3 4 5) -> (1 3 5)(0 2)
(0 1) -> (3 5)(2 4)(0 1)
)
Coset(
(0 1 2 3 4 5) -> (0 1 2 3 4 5)
(0 1) -> (0 1)
)
```
with the first element being the lone non-trivial element. The coset representative displayed here sends transpositions to a product of three transpositions, which is another well known property of these automorphisms. An automorphism representing this non-trivial element of $\mathrm{Out}(S_6)$ can be called via:
```python
exotic_auto = OutG[0].g
```
the 0th element is the coset, and the `.g` attribute selects the left factor of a `Coset` representing $gH$. There are of course 720 total automorphisms of $S_6$ which are not inner, and these can be called via `OutG[0][i]` where `i` ranges from 0 to 719. Now, `exotic_auto` can be called just like any other function to act on elements of $S_6$:
```python
exotic_auto(G[3])
```
which returns `(1 3 4)(0 5 2)`.