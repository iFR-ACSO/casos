# CaΣoS, the nonlinear sum-of-squares optimization suite

CaΣoS provides tools for symbolic polynomial expressions as well as parametrized convex and nonconvex sum-of-squares optimization problems, making use of the [CasADi](https://web.casadi.org) software for symbolic expressions, automatic differentiation, and numerical optimization.

### Install

The following requirements need to be met in order to use all functionalities of CaΣoS:
1. Download [CasADi v3.6.x](https://web.casadi.org/get/#body-36) and add it to your Matlab path.
2. Download and install at least one solver for semidefinite (conic) programs, and add the solver(s) to the Matlab Path.

   - Currently supported are [Mosek](https://www.mosek.com/downloads/) (v10.1), [SCS](https://www.cvxgrp.org/scs/install/matlab.html#matlab-install) (v3.2.4), and [SeDuMi](https://sedumi.ie.lehigh.edu/?page_id=58) (v1.3).
   - **Note:** Mosek requires a separate license. [Academic Licenses](https://www.mosek.com/products/academic-licenses/) are available.

3. Add the CaΣoS root folder (the one that contains the directory `+casos`) to your Matlab path.

> [!IMPORTANT]
> CaΣoS requires CasADi version 3.6.0 or newer.

## Polynomial expressions

The classes `casos.PD` and `casos.PS` implement polynomials of which the coefficients are constant doubles or can be symbolic expressions, respectively.

#### Polynomials of degree zero
Polynomials of degree zero correspond to constant or symbolic expressions without indeterminate variables.

```
casos.PD(M)
casos.PS(M)
```
creates a zero-degree polynomial which corresponds to the double or `casadi.DM`, or `casadi.SX` matrix `M`.

The following syntaxes for constant matrices are also supported by `casos.PS`:
```
casos.PD(m,n)
casos.PD.zeros(m,n)
casos.PD.zeros(n)
```
creates a zero-degree polynomial which corresponds to a `m × n` matrix (resp., a column vector with length `n`) of zeros.

```
casos.PD.ones(m,n)
casos.PD.ones(n)
```
creates a zero-degree polynomial which corresponds to a `m × n` matrix (resp., a column vector with length `n`) of ones.

```
casos.PD.eye(n)
```
creates a zero-degree polynomial which corresponds to the `n × n` identity matrix.

#### Indeterminate variables
CaΣoS distinguishes between *indeterminate* variables (symbols in a polynomial sense) and *symbolic* variables (variables in an optimization sense).

```
casos.Indeterminates('x',n)
casos.Indeterminates('x','y',...)
```
creates a tupel of `n` indeterminate variables; in the second case, `n` corresponds to the number of arguments.

Tupels of indeterminate variables can be converted into polynomials that correspond to vectors of indeterminate variables, and vice-versa. Moreover, indeterminate variables can be used in algebraic expressions to define polynomials with constant or symbolic coefficients, e.g.,
```
f = [-x(2); x(1) + (x(1)^2 - 1)*x(2)]
u = K*x
```
if `x` is a tuple of indeterminate variables and `K` is a double, `casadi.DM`, or `casadi.SX` matrix of suitable dimensions.

#### Monomial patterns
Monomial patterns describe the monomial terms a polynomial expression has.

```
monomials(x,deg)
```
creates a pattern of the `l` monomials in `x` with degree(s) in `deg`, where `x` must be a tuple of indeterminate variables (or a corresponding, vector-valed polynomial) and `deg` is a list of nonnegative integers; where `l` is the total number of such monomials.

Monomial patterns are defined by the class `casos.Sparsity`; if multiple patterns are concatenated to a vector or matrix, we call that a monomial *sparsity* pattern. The previous syntax is equivalent to
```
casos.Sparsity.scalar(x,deg)
```
while the following syntaxes can be used to create more complex, matrix-valued monomial sparsity patterns:
```
casos.Sparsity.dense(...,w)
casos.Sparsity.diag(...,w)
casos.Sparsity.band(...,w)
casos.Sparsity.banded(...,w)
casos.Sparsity.nonzeros(...,w)
casos.Sparsity.triplet(...,w)
```
where `...` denotes arguments to the equivalent function of `casadi.Sparsity` and `w` is either a scalar monomial pattern or a monomial sparsity pattern with length equal to the nonzero matrix entries.

Both (scalar) monomial patterns and monomial sparsity patterns can be used to define symbolic polynomials.

#### Polynomials with symbolic coefficients
The class `casos.PS` implements polynomials of which the coefficients can be symbolic expressions.
Unlike indeterminate variables, symbolic polynomials (variables) can be decision variables of an optimization problem. The following syntax creates polynomials which have symbolic variables as coefficients. In all of the following syntaxes, the first argument corresponds to the display name (resp., prefix) for the symbolic coeffcients. See [Casadi's `SX` symbolics](https://web.casadi.org/docs/#the-sx-symbolics) for details.

```
casos.PS.sym('c',w)
```
creates a scalar polynomial with symbolic coefficients and monomials in `w`, if `w` is a scalar monomial pattern; *or* creates a vector/matrix of polynomials with symbolic coefficients, size equal to the size of `w`, and its `(i,j)`-th entry having the monomial terms of `w(i,j)`.

```
casos.PS.sym('c',w,[m n])
casos.PS.sym('c',w,n)
```
creates a `m × n` matrix or a `n × 1` vector of polynomials with symbolic coefficients and each entry having the monomial terms in `w`, where `w` must be a scalar monomial pattern.

```
casos.PS.sym('c',[m n])
casos.PS.sym('c',n)
```
creates a `m × n` matrix or a `n × 1` vector of polynomials with degree zero; essentially, this is a symbolic matrix similar to `casadi.SX`.

```
casos.PS.sym(...,'gram')
```
where `...` denotes any of the syntaxes above, creates a scalar or matrix polynomial in Gram form, that is, with entries `p = z'*Q*z`, where `z` is the vector of monomials in `w` and `Q` is a quadratic symbolic matrix.

**Note:** We say that a polynomial is *symbolic* if and only if all of its (nonzero) coefficients are symbols in the sense of Casadi. All polynomials created with the syntaxes above, *including the Gram form,* are symbolic but the result of the notation 
```
casos.PS.sym('c',[1 2])*[x;x]
```
with `x = casos.PS('x')` would only be a symbolic *expression*. The queries `is_symbolic` and `is_symexpr` check whether a polynomial is a symbolic polynomial or a symbolic expression, respectively.

## Functions between polynomials

The class `casos.Function` provides functions of which the input and/or output are polynomial expressions.

```
casos.Function('f',{p1 ... pN},{q1 ... qM})
```
creates a function named `'f'` mapping the `M` outputs to `N` inputs. Outputs may be any expressions of types `casadi.DM`, `casadi.SX`, `casos.PD`, or `casos.PS` whereas inputs must be *symbolic* expressions of types `casadi.SX` or `casos.PS`.

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

## Sum-of-squares optimization

A polynomial sum-of-squares optimization problem takes the form

```math
\begin{array}{l c r}
  \min & F(\xi,\pi), & \xi = (\xi_\mathrm{l}, \xi_\mathrm{c}) \\
  \text{s.t.} & \gamma_\mathrm{lb} \unlhd G_\mathrm{l}(\xi,\pi) \unlhd \gamma_\mathrm{ub}, & G_\mathrm{c}(\xi,\pi) \in \mathcal K_c \\
  \text{and} & \xi_\mathrm{lb} \unlhd \xi_\mathrm{l} \unlhd \xi_\mathrm{ub}, & \xi_\mathrm{c} \in \mathcal K_x
\end{array}
```
where $F$ is a scalar-valued function, the constraints $G_\mathrm{l}$ and $G_\mathrm{c}$ take polynomial values, and $\unlhd$ denotes a coefficient-wise inequality; $\mathcal K_c$ and $\mathcal K_x$ are convex cones in the space of polynomials. The pairs of polynomials $(\xi_\mathrm{lb}, \gamma_\mathrm{lb})$ denote *lower bounds* and $(\xi_\mathrm{ub}, \gamma_\mathrm{ub})$ denote *upper bounds*.

#### Affine interface

A sum-of-squares problem is affine if $F$ is a linear (or quadratic) form in $\xi$ and $G = (G_\mathrm{l}, G_\mathrm{c})$ are affine functions in $\xi$.

```
S = casos.sossol('S','solver',struct('x',xi,'f',F,'g',G,'p',pi),opts)
```
initializes the SOS solver named `'S'` by relaxation to a convex optimization problem using the convex solver `'solver'`. See [Convex optimization](#convex-optimization) for supported solvers. Options are provided as structure `opts` including optional fields `opts.Kx` and `opts.Kc` describing, respectively, the cones $\mathcal K_x$ and $\mathcal K_c$. See [Polynomial cones](#polynomial-cones) for details.

```
sol = S('lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg)
```
evaluates the SOS solver `S` providing (optional) arguments to describe $\xi_\mathrm{lb}$, $\xi_\mathrm{ub}$ and $\gamma_\mathrm{lb}$, $\gamma_\mathrm{ub}$.

#### Quasiconvex interface

A quasiconvex sum-of-squares problems takes the form

```math
\begin{array}{l c r}
  \min & \pm t, & \xi = (\xi_\mathrm{l}, \xi_\mathrm{c}) \\
  \text{s.t.} & \gamma_\mathrm{lb} \unlhd G_\mathrm{l}(t,\xi,\pi) \unlhd \gamma_\mathrm{ub}, & G_\mathrm{c}(t,\xi,\pi) \in \mathcal K_c \\
  \text{and} & \xi_\mathrm{lb} \unlhd \xi_\mathrm{l} \unlhd \xi_\mathrm{ub}, & \xi_\mathrm{c} \in \mathcal K_x
\end{array}
```
where $t$ enters affinely into $G = (G_\mathrm{l}, G_\mathrm{c})$ and $G(t, \xi, \pi)$ is affine in $\xi$ for any $t \in \mathbb R$. For details on quasiconvex sum-of-squares programming, refer to [[SB2010]](#references).

```
S = casos.qcsossol('S','bisection',struct('x',xi,'f',±t,'g',G,'p',pi),opts)
```
initializes the quasiconvex SOS solver named `'S'` by bisection over convex sum-of-squares optimization problems. The options structure `opts` includes optional fields `opts.Kx` and `opts.Kc` describing $\mathcal K_x$ and $\mathcal K_c$, respectively. See [Polynomial cones](#polynomial-cones) for details.

#### Polynomial cones

The options `K_` define the (convex) polynomial cones $\mathcal K$ as well as the number of coefficient-wise constraints in the SOS problems. Each option takes a structure `K` as value which can have the following fields:

- `K.lin` : number of coefficient-wise constraints; corresponds to the dimension of $G_\mathrm{l}$ or $\xi_\mathrm{l} \in \mathbb R[x]^l$;
- `K.sos` : number of sum-of-squares constraints, that is, $\mathcal K = \Sigma[x]^s$; corresponds to the dimension of $g_\mathrm{c}$ or $\xi_\mathrm{c}$;
- `K.dsos` : number of diagonally dominant sum-of-squares constraints, that is, $\mathcal K = DSOS[x]^s$;
- `K.sdsos` : number of scaled diagonally dominant sum-of-squares constraints, that is, $\mathcal K = SDSOS[x]^s$;
- no further cones are currently supported;

by default (if the option `K_` is omitted), only coefficient-wise constraints are enforced.

### Convex optimization

To solve polynomial sum-of-squares optimization problems via semidefinite programming, CaΣoS provides both a high-level and low-level interface for convex optimization. These interfaces extend CasADi's [`qpsol`](https://web.casadi.org/docs/#high-level-interface) and [`conic`](https://web.casadi.org/docs/#low-level-interface) syntax. The convex solvers SeDuMi, MOSEK,[^2] and SCS are supported (but, unlike CasADi, the solvers must be installed separately and be accessible through MATLAB). 

> [!NOTE]
> SeDuMi does not support quadratic cost functions.

[^2]: Using the MOSEK interface requires a separate licence; please see [MOSEK's homepage](https://www.mosek.com/resources/getting-started/) for details.

#### High-level interface

The high-level interface solves convex problems of the form

```math
\begin{array}{l c r}
  \min & f(x,p), & x = (x_\mathrm{l}, x_\mathrm{c}) \\
  \text{s.t.} & g_\mathrm{lb} \leq g_\mathrm{l}(x,p) \leq g_\mathrm{ub}, & g_\mathrm{c}(x,p) \succeq_{\mathcal K_c} g_\mathrm{cb} \\
  \text{and} & x_\mathrm{lb} \leq x_\mathrm{l} \leq x_\mathrm{ub}, & x_\mathrm{c} \succeq_{\mathcal{K}_x} x_\mathrm{cb}
\end{array}
```
where $f$ is a convex quadratic function in $x$, the constraints $g_\mathrm{l}$ and $g_\mathrm{c}$ are affine in $x$, and $\succeq_\mathcal{K}$ denotes the order induced by the convex cone $\mathcal K$. Moreover, the pairs $(x_\mathrm{lb}, g_\mathrm{lb})$ denote *lower bounds*, $(x_\mathrm{ub}, g_\mathrm{ub})$ denote *upper bounds*, and $(x_\mathrm{cb}, g_\mathrm{cb})$ denote *conic bounds*.

```
S = casos.sdpsol('S','solver',struct('x',x,'f',f,'g',g,'p',p),opts)
```
initializes the SDP solver named `'S'` using the convex solver `'solver'`. Options are provided as structure `opts` including optional fields `opts.Kx` and `opts.Kc` describing, respectively, the cones $\mathcal K_x$ and $\mathcal K_c$. See [Convex cones](#convex-cones) for details.

```
sol = S('lbx',lbx,'ubx',ubx,'cbx',cbx,'lbg',lbg,'ubg',ubg,'cbg',cbg)
```
evaluates the SDP solver `S` providing (optional) arguments to describe $x_\mathrm{lb}$, $x_\mathrm{ub}$, $x_\mathrm{cb}$ and $g_\mathrm{lb}$, $g_\mathrm{ub}$, $g_\mathrm{cb}$.

#### Low-level interface

The low-level interface solves conic problems of the form

```math
\begin{array}{l c r}
  \min & \frac{1}{2} x^\top H x + g^\top x, & x = (x_\mathrm{l}, x_\mathrm{c}) \\
  \text{s.t.} & a_\mathrm{lb} \leq A_\mathrm{l} \, x \leq a_\mathrm{ub}, & A_\mathrm{c} \, x \succeq_{\mathcal K_c} a_\mathrm{cb} \\
  \text{and} & x_\mathrm{lb} \leq x_\mathrm{l} \leq x_\mathrm{ub}, & x_\mathrm{c} \succeq_{\mathcal{K}_x} x_\mathrm{cb}
\end{array}
```
where $\succeq_\mathcal{K}$ denotes the order induced by the convex cone $\mathcal K$. Moreover, the pairs $(x_\mathrm{lb}, a_\mathrm{lb})$ denote *lower bounds*, $(x_\mathrm{ub}, a_\mathrm{ub})$ denote *upper bounds*, and $(x_\mathrm{cb}, a_\mathrm{cb})$ denote *conic bounds*.

```
S = casos.conic('S','solver',struct('h',hs,'a',as),opts)
```
initializes the conic solver named `'S'` using the convex solver `'solver'`, where `hs` and `as` are sparsity patterns for $H$ and $A = (A_\mathrm{l}, A_{\mathrm c})$. Options are provided as structure `opts` including optional fields `opts.Kx` and `opts.Kc` describing, respectively, the cones $\mathcal K_x$ and $\mathcal K_c$. See [Convex cones](#convex-cones) for details.

```
sol = S('h',h,'g',g,'a',a,'lba',lba,'uba',uba,'cba',cba,'lbx',lbx,'ubx',ubx,'cbx',cbx)
```
evaluates the conic solver `S` providing (optional) arguments to describe $H$, $A$, and $g$ as well as $a_\mathrm{lb}$, $a_\mathrm{ub}$, $a_\mathrm{cb}$ and $x_\mathrm{lb}$, $x_\mathrm{ub}$, $x_\mathrm{cb}$.

#### Convex cones

The options `K_` define the convex cones $\mathcal K$ as well as the number of linear constraints in SDP or conic problems. Each option takes a structure `K` as value which can have the following fields:

- `K.lin` : number of linear constraints; corresponds to the first dimension of $A_\mathrm{l}$, the dimension of $g_\mathrm{l}$, or the dimension of $x_\mathrm{l}$;
- `K.psd` : list $(s_1, \ldots, s_k)$ of semidefinite cone dimensions, that is, $\mathcal K = \mathbb S_{s_1}^{+} \times \cdots \times \mathbb S_{s_k}^{+}$; the total number of SDP cone constraints is equal to $\sum_i s_i^2$;
- `K.dd` : list $(d_1, \ldots, d_k)$ of diagonally dominant cone dimensions, that is, $\mathcal K =  DD_{d_1} \times \cdots \times DD_{d_k}$ with $DD_{d}$ being the set of $d \times d$ diagonally dominant matrices; the total number of DD cone constraints is equal to $\sum_i d_i^2$;
- `K.sdd` : list $(d_1, \ldots, d_k)$ of scaled diagonally dominant cone dimensions, that is, $\mathcal K =  SDD_{s_1} \times \cdots \times SDD_{s_k}$ with $SDD_{s}$ being the set of $s \times s$ scaled diagonally dominant matrices; the total number of SDD cone constraints is equal to $\sum_i d_i^2$;
- further cones, e.g., the Lorentz (or second-order) cone, are supported depending on the convex solver; the total number of *all* cone constraints corresponds to the first dimension of $A_\mathrm{c}$, the dimension of $g_\mathrm{c}$, or the dimension of $x_\mathrm{c}$;

by default (if the option `K_` is omitted), only linear constraints are enforced.

## Transitioning

The following comparison is supposed to ease the transition from other sum-of-squares or general optimization toolboxes to CaΣoS.

> [!NOTE]
> This section only shows a subset of the CaΣoS interface. For full details, please see the descriptions above.

#### SOSOPT

| SOSOPT | Description | CaΣoS |
|--------|-------------|-------|
| `polynomial(0)` | Constant polynomial. | `casos.PS(0)` |
| `pvar('x')` | Scalar indeterminate variable. | `casos.PS('x')` |
| `pvar('q')` | Scalar decision variable. | `casos.PS.sym('q')` |
| `mpvar('x',n,m)` | Matrix of indeterminate variables. | `casos.PS('x',n,m)` |
| `mpvar('Q',n,m)` | Matrix decision variable. | `casos.PS.sym('Q',n,m)` |
| `monomials(x,deg)` | Vector of monomials. | `monomials(x,deg)` |
| `polydecvar('c',z)` | Polynomial decision variable $c^\top z$. | `casos.PS.sym('c',z)` |
| `sosdecvar('Q',z)` | Gram decision variable $z^\top Q z$. | `casos.PS.sym('Q',z,'gram')` |
| `jacobian(f,x)` | Partial derivative w.r.t. indeterminates. | `nabla(f,x)` |
| `jacobian(p,q)` | Partial derivative w.r.t. symbolic variables. | *Not yet supported* |
| `constr = (expr >= 0)` | Sum-of-squares expression constraint. | `sos.g = expr; opts.Kc.sos = 1` |
| `constr = (svar >= 0)` | Sum-of-squares variable constraint <br/> (requires Gram variable). | `sos.x = svar; opts.Kx.sos = 1` |
| `constr = (p == q)` | Polynomial expression equality. | `sos.g = (p - q); opts.Kc.lin = 1` <br/> `lbg = 0; ubg = 0` |
| `constr = (q <= 1)` | Scalar variable inequality. | `sos.x = q; opts.Kx.lin = 1` <br/> `lbx = -inf; ubx = 1` |
| `sosopt(constr,x)` | Sum-of-squares feasibility. | `S = casos.sossol('S','solver',sos,opts)` |
| `sosopt(constr,x,obj)` | Sum-of-squares optimization. | `sos.f = obj` <br/> `S = casos.sossol('S','solver',sos,opts)` |
| `[info,dopt] = sosopt(...)` | Solve affine problem. | `S = casos.sossol('S','solver',sos,opts)` <br/> `sol = S(...)` |
| `info.feas` | Retrieve feasibility info. | `S.stats.UNIFIED_RETUR_STATUS` |
| `info.obj` | Retrieve optimal value. | `sol.f` |
| `gsosopt(constr,x,obj)` | Quasi-convex optimization (bisection). | `sos.f = obj` <br/> `S = casos.qcsossol('S','bisection',sos,opts)` |
| `[info,dopt] = gsosopt(...)` | Solve quasi-convex problem. | `sol = S(...)` |
| `info.tbnds(2)` | Retrieve upper bound on optimal value. | `sol.f` |
| `subs(s,dopt)` | Retrieve optimal solution (variable). | `sol.x` |
| `subs(p,dopt)` | Retrieve optimal solution (expression). | `sol.g` |

## References

[SB2010]: Seiler, P., Balas, G.J.: Quasiconvex sum-of-squares programming. *49th IEEE Conference on Decision and Control*, pp. 3337–3342, Atlanta, GA (2010). [10.1109/CDC.2010.5717672](https://doi.org/10.1109/CDC.2010.5717672).
