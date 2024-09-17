# Multipoly: A Toolbox for Multivariate Polynomials

Version 2.00, 23 November 2010

## Creating polynomial objects
   [`PVAR`](pvar.m)         - Construct a polynomial variable
   
   [`MPVAR`](mpvar.m)        - Construct a matrix or vector polynomial variable
   
   [`POLYNOMIAL`](@polynomial/polynomial.m)   - Construct a polynomial object
   
   [`MONOMIALS`](monomials.m)    - Construct list of monomials
   
   [`PDATAFIT`](@polynomial/pdatafit.m)     - Compute a polynomial least squares fit to data
   
   [`PFUNCTIONFIT`](@polynomial/pfunctionfit.m) - Compute a polynomial least squares fit to a function
   
## Simulink:
   [POLYLIB.MDL](polylib.mdl)  - Simulink block for polynomial objects

## Polynomial plotting:
   [`PCONTOUR`](@polynomial/pcontour.m)     - Plot 2d polynomial contours
   
   [`PCONTOUR3`](@polynomial/pcontour3.m)    - Plot 3d polynomial contours

## Polynomial functions:
   [`POLY2BASIS`](@polynomial/poly2basis.m)   - Project polynomial onto a basis of monomials
   
   [`PLINEARIZE`](@polynomial/plinearize.m)   - Linearize a vector polynomial function
   
   [`PTRIM`](@polynomial/ptrim.m)        - Find trim conditions (equilibria) for a polynomial dynamical system
   
   [`PVOLUME`](@polynomial/pvolume.m)      - Estimate the volume of a polynomial set
   
   [`PSAMPLE`](@polynomial/psample.m)      - Draw random samples from a polynomial set
   
   [`PSIM`](@polynomial/psim.m)         - Simulate a polynomial dynamical system
   
   [`PPLANESIM`](@polynomial/pplanesim.m)    - Plot the phase plane for a polynomial dynamical system
   
   [`INT`](@polynomial/int.m)          - Element-by-element integration of a polynomial
   
   [`DIFF`](@polynomial/diff.m)         - Element-by-element differentiation of a polynomial
   
   [`JACOBIAN`](@polynomial/jacobian.m)     - Compute Jacobian matrix of a polynomial vector
   
   [`COLLECT`](@polynomial/collect.m)      - Collect coefficients of specified variables
   
   [`SUBS`](@polynomial/subs.m)         - Symbolic substitution
   
   [`CLEANPOLY`](@polynomial/cleanpoly.m)    - Remove terms based on value of coefficient and degree

## Polynomial characteristics:
   [`ISDOUBLE`](isdouble.m)     - True for arrays of doubles
   
   [`ISPVAR`](ispvar.m)       - True for arrays of pvars
   
   [`ISMONOM`](ismonom.m)      - True for arrays of monomials
   
   [`ISEMPTY`](@polynomial/isempty.m)      - True for empty monomials
   
   [`ISEQUAL`](@polynomial/isequal.m)      - Element by element polynomial comparisons
   
   [`SIZE`](@polynomial/size.m)         - Size of a polynomial matrix
   
   [`LENGTH`](@polynomial/length.m)       - Length of a polynomial matrix
   
   [`FIELDNAMES`](@polynomial/fieldnames.m)   - Get properties of a polynomial object

## Conversions:
   [`P2S`](p2s.m)          - Convert from multipoly to symbolic toolbox
   
   [`S2P`](s2p.m)          - Convert from symbolic toolbox to multipoly
   
   [`DOUBLE`](@polynomial/double.m)       - Convert constant polynomial to a double
   
   [`CHAR`](@polynomial/char.m)         - Converts a polynomial to its string representation.

## Overloaded arithmetic operations:
   [`PLUS, +`](@polynomial/plus.m)      - Add polynomials
   
   [`MINUS, -`](@polynomial/minus.m)     - Subtract polynomials
   
   [`MTIMES, *`](@polynomial/mtimes.m)    - Multiply polynomials
   
   [`MPOWER, ^`](@polynomial/mpower.m)    - Power of a polynomial
   
   [`HORZCAT, [,]`](@polynomial/horzcat.m) - Horizontal concatentation of polynomials
   
   [`VERTCAT, [;]`](@polynomial/vertcat.m) - Vertical concatentation of polynomials
   
   [`DIAG`](@polynomial/diag.m)         - Diagonal poly matrices and diagonals of poly matrices
   
   [`TRIL`](@polynomial/tril.m)         - Extract lower triangular part of a polynomial matrix
   
   [`TRIU`](@polynomial/triu.m)         - Extract upper triangular part of a polynomial matrix
   
   [`BLKDIAG`](@polynomial/blkdiag.m)      - Block diagonal concatenation of polynomial matrices
   
   [`CTRANSPOSE, '`](@polynomial/ctranspose.m) - Complex-conjugate transpose of a polynomial
   
   [`TRANSPOSE, .'`](@polynomial/transpose.m) - Non-conjugate transpose of a polynomial
   
   [`RESHAPE`](@polynomial/reshape.m)      - Reshape a polynomial matrix
   
   [`REPMAT`](@polynomial/repmat.m)       - Replicate and tile an array of polynomials
   
   [`UPLUS`](@polynomial/uplus.m)        - Unary plus of a polynomial
   
   [`UMINUS`](@polynomial/uminus.m)       - Unary minus of a polynomial
   
   [`TIMES, .*`](@polynomial/times.m)    - Element-by-element multiply of polynomials
   
   [`POWER, .^`](@polynomial/power.m)    - Element-by-element power of a polynomial
   
   [`SUM`](@polynomial/sum.m)          - Sum of the elements of a polynomial array
   
   [`PROD`](@polynomial/prod.m)         - Product of the elements of a polynomial array
   
   [`TRACE`](@polynomial/trace.m)        - Sum of the diagonal elements
   
   [`DET`](@polynomial/det.m)          - Determinant of a polynomial matrix

## More information
See [multipolydoc.pdf](doc/multipolydoc.pdf) for a full documentation.
