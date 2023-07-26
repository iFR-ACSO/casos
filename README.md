# CaΣoS: CasADi-based sum-of-squares optimization suite[^1]

## Polynomial expressions

The class `casos.PS` implements polynomials of which the coefficients can be symbolic expressions.

#### Polynomials of degree zero

```
casos.PS(M)
```
creates a zero-degree polynomial which corresponds to the double, `casadi.DM`, or `casadi.SX` matrix `M`.

```
casos.PS(m,n)
casos.PS.zeros(m,n)
casos.PS.zeros(n)
```
creates a zero-degree polynomial which corresponds to a `m`x`n` matrix (resp., a square matrix with length `n`) of zeros.

```
casos.PS.ones(m,n)
casos.PS.ones(n)
```
creates a zero-degree polynomial which corresponds to a `m`x`n` matrix (resp., a square matrix with length `n`) of ones.

```
casos.PS.eye(n)
```
creates a zero-degree polynomial which corresponds to the `n`x`n` identity matrix.

#### Indeterminate variables & monomials

```
casos.PS('x','y',...)
```
creates a `n`x1 vector of indeterminate variables, where `n` corresponds to the number of arguments.

```
casos.PS('x',m,n)
casos.PS('x',n)
```
creates a `m`x`n` matrix (resp., a square matrix) of indeterminate variables.

```
monomials(x,deg)
```
creates a `l`x1 vector of all monomials in `x` with degree in `deg`, where `x` must be a vector of indeterminate variables and `deg` is vector of nonnegative integers; where `l` is the total number of such monomials.

```
monomials(p)
```
creates a `l`x1 vector of all monomials in the polynomial `p`; where `l` is the total number of monomials in `p`.

#### Polynomials with symbolic coefficients

```
casos.PS.sym('c',w)
```
creates a scalar polynomial with symbolic coefficients and monomials in `w`, where `w` must be a vector of monomials.

```
casos.PS.sym('c',w,[m n])
casos.PS.sym('c',w,n)
```
creates a `m`x`n` matrix (resp., a square matrix with length `n`) of polynomials with symbolic coefficients and monomials in `w`, where `w` must be a vector of monomials.

## Functions between polynomials

The class `casos.Function` provides functions of which the input and/or output are polynomial expressions.

```
casos.Function('f',{p1 ... pN},{q1 ... qM})
```
creates a function named `'f'` mapping the `M` outputs to `N` inputs. Outputs may be any expressions of types `casadi.DM`, `casadi.SX`, or `casos.PS` whereas inputs must be *symbolic* expressions of types `casadi.SX` or `casos.PS`.

```
casos.Function('f',{p1 ... pN},{q1 ... qM},{'a1' ... 'aN'},{'b1' ... 'bM'})
```
creates a function named `'f'` as described above but also assigns names to input and outputs.

#### Calling functions

Suppose `f` is a `casos.Function` object.

```
[s1,...,sM] = f(r1,...,rN)
```
calls the function `f` with arguments `r1` through `rN` assigned to its inputs and returns the output values `s1` through `sM`. 

If the function `f` has named inputs and outputs, the following call syntax is also possible.

```
out = f('a1',r1,...,'aN',rN)
```
calls the function `f` with arguments `r1` through `rN` assigned to inputs with names `'a1'` through `'aN'`, respectively. In this case, the functionn call returns a structure with fields `out.b1` through `out.bM` with corresponding output values.


[^1]: CaΣoS has been neither supported nor endorsed by CasADi or any of its affilitiates.
