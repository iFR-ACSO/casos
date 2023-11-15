# CaΣoS: _CasADi-based sum-of-squares optimization suite_[^1]

CaΣoS provides a symbolic framework for convex and nonconvex sum-of-squares problems, making use of the [CasADi](https://web.casadi.org) software for symbolic expressions, automatic differentiation, and numerical optimization.

## Polynomial expressions

The class `casos.PS` implements polynomials of which the coefficients can be symbolic expressions.

#### Polynomials of degree zero
Polynomials of degree zero correspond to constant or symbolic expressions without indeterminate variables.

```
casos.PS(M)
```
creates a zero-degree polynomial which corresponds to the double, `casadi.DM`, or `casadi.SX` matrix `M`.

```
casos.PS(m,n)
casos.PS.zeros(m,n)
casos.PS.zeros(n)
```
creates a zero-degree polynomial which corresponds to a `m × n` matrix (resp., a column vector with length `n`) of zeros.

```
casos.PS.ones(m,n)
casos.PS.ones(n)
```
creates a zero-degree polynomial which corresponds to a `m × n` matrix (resp., a column vector with length `n`) of ones.

```
casos.PS.eye(n)
```
creates a zero-degree polynomial which corresponds to the `n × n` identity matrix.

#### Indeterminate variables & monomials
CaΣoS distinguishes between *indeterminate* variables (symbols in a polynomial sense) and *symbolic* variables (variables in an optimization sense). The following syntax creates polynomials that correspond to (vectors of) indeterminate degrees and monomial expressions.

```
casos.PS('x','y',...)
```
creates a `n × 1` vector of indeterminate variables, where `n` corresponds to the number of arguments.

```
casos.PS('x',m,n)
casos.PS('x',n)
```
creates a `m × n` matrix (resp., a square matrix) of indeterminate variables.

```
monomials(x,deg)
```
creates a `l × 1` vector of all monomials in `x` with degree(s) in `deg`, where `x` must be a vector of indeterminate variables and `deg` is vector of nonnegative integers; where `l` is the total number of such monomials.

```
monomials(p)
```
creates a `l × 1` vector of all monomials in the polynomial `p`; where `l` is the total number of monomials in `p`.

#### Polynomials with symbolic coefficients
Unlike indeterminate variables, symbolic variables can be decision variables of an optimization problem. The following syntax creates polynomials which have symbolic variables as coefficients. In all of the following syntaxes, the first argument corresponds to the display name (resp., prefix) for the symbolic coeffcients. See [Casadi's `SX` symbolics](https://web.casadi.org/docs/#the-sx-symbolics) for details.

```
casos.PS.sym('c',w)
```
creates a scalar polynomial with symbolic coefficients and monomials in `w`, where `w` must be a vector of monomials.

```
casos.PS.sym('c',w,[m n])
casos.PS.sym('c',w,n)
```
creates a `m × n` matrix (resp., a square matrix with length `n`) of polynomials with symbolic coefficients and monomials in `w`, where `w` must be a vector of monomials.

```
casos.PS.sym('c',[m n])
casos.PS.sym('c',n)
```
creates a `m × n` matrix (resp., a square matrix with length `n`) of polynomials of degree zero; essentially, this is a symbolic matrix similar to `casadi.SX`.

```
casos.PS.sym(...,'gram')
```
where `...` denotes any of the syntaxes above, creates a scalar or matrix polynomial in Gram form, that is, with entries `p = z'*Q*z`, where `z` is the vector of monomials in `w` and `Q` is a quadratic symbolic matrix.

**Note 1:** The syntax `casos.PS.sym('c',w,...)` creates polynomials with monomials *in* `w` and is therefore equivalent to `casos.PS.sym('c',monomials(w),...)`. In consequence, all of the following syntaxes yield the same result:
```
casos.PS.sym('c',[1;w])
casos.PS.sym('c',[w;1])
casos.PS.sym('c',[1;w;w])
```

**Note 2:** We say that a polynomial is *symbolic* if and only if all of its (nonzero) coefficients are symbols in the sense of Casadi. Except for the Gram form, all polynomials created with the syntaxes above are symbolic but the result of the notation 
```
casos.PS.sym('c',[1 2])*[x;x]
```
with `x = casos.PS('x')` would only be a symbolic *expression*. The same is true for the Gram form syntax because of the symmetries in the Gram matrix expression. The queries `is_symbolic`, `is_symexpr`, and `is_symgram` check whether a polynomial is a symbolic polynomial, a symbolic expression, or a symbolic Gram form, respectively.

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
