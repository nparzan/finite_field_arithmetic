# Finite Field Arithmetic

A simple Python module for performing finite field arithmetic.

![GIF demo](media/demo.gif)

## Usage
This module implements three algebraic elements:
* Fields
* Polynomials
* Field elements

In order to load the module simple import it:
```python
from finite_field_arith import *
```

### Fields

The general syntax for creating a field is:

```python
F = Field(p, dim, irr)
```

This is the construction of the field of size `p^dim` where `p` is some prime and `dim` is a positive integer. If `dim > 1` then `irr` is an irreducible polynomial in `F_p[x]` of degree `n`.

For example, to construct `F_2`:
```python
F2 = Field(2, 1, None)
```

Given an irreducible polynomial of degree `3` over `F_2`, say `irr = x^3 + x + 1` one can construct `F_8` simply by calling:

```python
F8 = Field(2, 3, irr)
```
### Polynomials
A polynomial `f` of degree `k` is determined by:
* Its underlying field `F`
* Its list of coefficients `alpha_0, ... , alpha_k`, all elements of `F`

Given a field `F` and a list of coefficients `coefs` from list we can construct the polynomial:
```python
f = Polynomial(coefs, F)
```
### Field Elements
A field element `alpha` in `F = G[x] % irr` is determined by the polynomial representation `f_alpha` and its respective field `F`. Note that while `f_alpha` is an element in `G[x]` (i.e., the coefficient are element from `G`) the field element `alpha` is an element of `F` and therefore is defined as follows:
<!-- since the polynomial representation of `alpha` already encapsulates `F` the construction of `alpha` is dependent on the polynomial alone. -->
<!-- That is, if `f` is a polynomial in some field `F` then: -->
```python
f = Polynomial(coefs, G)
alpha = Element(f, F)
```
Constructs the appropriate field element `f % F.irr`

### Full example
The following snippet constructs `F_2, F_8 = F_2[x] % x^3 + x + 1`, the polynomial `f = 1 + x^2` over `F_2[x]` and the field element `alpha_f = 1 + x^2 % (x^3 + x + 1)`:
```python
F2 = Field(2, 1, None)
irr = Polynomial([1, 1, 0, 1], F2)
F8 = Field(2, 3, irr)
f = Polynomial([1, 0, 1], F2)
f_in_field = Element(f, F8)
```

## Functionality
In what follows, let `F` be a field, `f, g, h` polynomials and `a, b` field elements. Things you can do with this module include:

### Polynomial arithmetic:

```python
t = f * g #product
t = f + g #addition (and subtraction)
t = f / g #quotient 
t = f % g #remainder
```

### Field operations:
```python
a_gen = a.generated_subgroup() #set of elements generated by powers of a
g = Element.draw_generator(F, halt = k) #draw random element in F and return it if it generates the entire field. Give up after k attempts (by default halt is unbounded)
```

## Todo
The following is a list of features to be implemented in the future:
* Complete docstrings: very partial atm.
* Extension fields: in the current implementation a field of size `p^k` for `k > 1` must be constructed as a degree `k` extension of `F_p`. Future implementation should allow, for example, creating `G = F_p[x] % irr_1` and `H = G[y] % irr_2` where `irr_1` is in `F_p[x]` and `irr_2` is in `G[y]`.

## License
This project is distributed under the Apache license version 2.0 (see the LICENSE file in the project root).
