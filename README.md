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
CaΣoS distinguishes between *indeterminate* variables (symbols in a polynomial sense) and *symbolic* variables (variables in an optimization sense). The following syntax creates polynomials that correspond to (vectors of) indeterminate variables and monomial expressions.

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

## Convex optimization

To solve polynomial sum-of-squares optimization problems via semidefinite programming, CaΣoS provides both a high-level and low-level interface for convex optimization. These interfaces extend CasADi's [`qpsol`](https://web.casadi.org/docs/#high-level-interface) and [`conic`](https://web.casadi.org/docs/#low-level-interface) syntax. Currently, the convex solvers SeDuMi and SCS are supported (but, unlike CasADi, the solvers must be installed separately and be accessible through MATLAB). Note that SeDuMi does not support quadratic cost functions.

#### High-level interface

The high-level interface solves convex problems of the form

```math
\begin{array}{l c r}
  \min & f(x,p), & x = (x_\mathrm{l}, x_\mathrm{c}) \\
  \text{s.t.} & g_\mathrm{lb} \leq g_\mathrm{l}(x,p) \leq g_\mathrm{ub}, & g_\mathrm{c}(x,p) \succeq_{\mathcal{K}_g} g_\mathrm{cb} \\
  \text{and} & x_\mathrm{lb} \leq x_\mathrm{l} \leq x_\mathrm{ub}, & x_\mathrm{c} \succeq_{\mathcal{K}_x} x_\mathrm{cb}
\end{array}
```
where $f$ is a convex quadratic function in $x$, the constraints $g_\mathrm{l}$ and $g_\mathrm{c}$ are affine in $x$, and $\succeq_\mathcal{K}$ denotes the order induced by the convex cone $\mathcal K$. Moreover, the pairs $(x_\mathrm{lb}, g_\mathrm{lb})$ denote *lower bounds*, $(x_\mathrm{ub}, g_\mathrm{ub})$ denote *upper bounds*, and $(x_\mathrm{cb}, g_\mathrm{cb})$ denote *conic bounds*.

```
S = casos.sdpsol('S','solver',struct('x',x,'f',f,'g',g,'p',p),opts)
```
initializes the SDP solver named `'S'` using the convex solver `'solver'`. Options are provided as structure `opts` including optional fields `opts.Kx` and `opts.Kg` describing, respectively, the cones $\mathcal K_x$ and $\mathcal K_g$. See [Convex cones](#convex-cones) for details.

```
sol = S('lbx',lbx,'ubx',ubx,'cbx',cbx,'lbg',lbg,'ubg',ubg,'cbg',cbg)
```
evaluates the SDP solver `S` providing (optional) arguments to describe $x_\mathrm{lb}$, $x_\mathrm{ub}$, $x_\mathrm{cb}$ and $g_\mathrm{lb}$, $g_\mathrm{ub}$, $g_\mathrm{cb}$.

#### Low-level interface

The low-level interface solves conic problems of the form

```math
\begin{array}{l c r}
  \min & \frac{1}{2} x^\top H x + g^\top x, & x = (x_\mathrm{l}, x_\mathrm{c}) \\
  \text{s.t.} & a_\mathrm{lb} \leq A_\mathrm{l} \, x \leq a_\mathrm{ub}, & A_\mathrm{c} \, x \succeq_{\mathcal{K}_A} a_\mathrm{cb} \\
  \text{and} & x_\mathrm{lb} \leq x_\mathrm{l} \leq x_\mathrm{ub}, & x_\mathrm{c} \succeq_{\mathcal{K}_x} x_\mathrm{cb}
\end{array}
```
where $\succeq_\mathcal{K}$ denotes the order induced by the convex cone $\mathcal K$. Moreover, the pairs $(x_\mathrm{lb}, a_\mathrm{lb})$ denote *lower bounds*, $(x_\mathrm{ub}, a_\mathrm{ub})$ denote *upper bounds*, and $(x_\mathrm{cb}, a_\mathrm{cb})$ denote *conic bounds*.

```
S = casos.conic('S','solver',struct('h',hs,'a',as),opts)
```
initializes the conic solver named `'S'` using the convex solver `'solver'`, where `hs` and `as` are sparsity patterns for $H$ and $A = (A_\mathrm{l}, A_{\mathrm c})$. Options are provided as structure `opts` including optional fields `opts.Kx` and `opts.Ka` describing, respectively, the cones $\mathcal K_x$ and $\mathcal K_a$. See [Convex cones](#convex-cones) for details.

```
sol = S('h',h,'g',g,'a',a,'lba',lba,'uba',uba,'cba',cba,'lbx',lbx,'ubx',ubx,'cbx',cbx)
```
evaluates the conic solver `S` providing (optional) arguments to describe $H$, $A$, and $g$ as well as $a_\mathrm{lb}$, $a_\mathrm{ub}$, $a_\mathrm{cb}$ and $x_\mathrm{lb}$, $x_\mathrm{ub}$, $x_\mathrm{cb}$.

#### Convex cones

The options `K_` define the convex cones $\mathcal K$ as well as the number of linear constraints in SDP or conic problems. Each option takes a structure `K` as value which can have the following fields:

- `K.l` : number of linear constraints; corresponds to the first dimension of $A_\mathrm{l}$ or dimension of $g_\mathrm{l}$;
- `K.s` : vector $(s_1, \ldots, s_k)$ of semidefinite cone dimensions, that is, $\mathcal K = \mathbb S^+_{s_1} \times \cdots \times \mathbb S^+_{s_k}$; the total number of SDP cone constraints is equal to $\sum_i s_i^2$;
- further cones, e.g., the Lorentz (or second-order) cone, are supported depending on the convex solver; the total number of *all* cone constraints corresponds to the first dimension of $A_\mathrm{c}$ or dimension of $g_\mathrm{c}$;

by default (if the option `K_` is omitted), only linear constraints are enforced.

[^1]: CaΣoS has been neither supported nor endorsed by CasADi or any of its affilitiates.
