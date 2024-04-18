The PS class
============

Polynomials
-----------

One polynomial is mathematically defined as

:math:`p(x) = \sum_a \underbrace{c_a x^a}_{\textrm{monomial}} = \sum_a c_a x_1^{a_1}x_2^{a_2}\cdots x_n^{a_n}`

where:
 
- :math:`x_i = (x_1, \cdots, x_n)` represent the indeterminate variables.
- :math:`a = (a_1, \cdots, a_n)` are the indices/degrees.
- :math:`c = (c_1, \cdots, c_m)` are coefficients for each monomial.

**Example:**

Consider the polynomial 

:math:`p(x_1, x_2) = 3x_1^2 x_2 + 5x_1 x_2^3 - 2x_2^2`

This is a polynomial with

- :math:`x=(x_1, x_2)`
- :math:`c=(3, 5, -2)^{\top}` = :code:`PS.coeffs`
- :math:`a = \begin{bmatrix} 2 & 1 \\ 5 & 3 \\ 0 & 2 \end{bmatrix}` = :code:`PS.degmat`



Polynomial expressions
----------------------

The class ``casos.PS`` implements polynomials of which the coefficients can be symbolic expressions.

Polynomials of degree zero
~~~~~~~~~~~~~~~~~~~~~~~~~~

Polynomials of degree zero correspond to constant or symbolic expressions without indeterminate variables.

.. code-block:: matlab

    casos.PS(M)

creates a zero-degree polynomial which corresponds to the double, ``casadi.DM``, or ``casadi.SX`` matrix ``M``.

.. code-block:: matlab

    casos.PS(m,n)
    casos.PS.zeros(m,n)
    casos.PS.zeros(n)

creates a zero-degree polynomial which corresponds to a ``m × n`` matrix (resp., a column vector with length ``n``) of zeros.

.. code-block:: matlab

    casos.PS.ones(m,n)
    casos.PS.ones(n)

creates a zero-degree polynomial which corresponds to a ``m × n`` matrix (resp., a column vector with length ``n``) of ones.

.. code-block:: matlab

    casos.PS.eye(n)

creates a zero-degree polynomial which corresponds to the ``n × n`` identity matrix.

Indeterminate variables & monomials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CaΣoS distinguishes between *indeterminate* variables (symbols in a polynomial sense) and *symbolic* variables (variables in an optimization sense). The following syntax creates polynomials that correspond to (vectors of) indeterminate variables and monomial expressions.

.. code-block:: matlab

    casos.PS('x','y',...)

creates a ``n × 1`` vector of indeterminate variables, where ``n`` corresponds to the number of arguments.

.. code-block:: matlab

    casos.PS('x',m,n)
    casos.PS('x',n)

creates a ``m × n`` matrix (resp., a square matrix) of indeterminate variables.

.. code-block:: matlab

    monomials(x,deg)

creates a ``l × 1`` vector of all monomials in ``x`` with degree(s) in ``deg``, where ``x`` must be a vector of indeterminate variables and ``deg`` is vector of nonnegative integers; where ``l`` is the total number of such monomials.

.. code-block:: matlab

    monomials(p)

creates a ``l × 1`` vector of all monomials in the polynomial ``p``; where ``l`` is the total number of monomials in ``p``.

Polynomials with symbolic coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unlike indeterminate variables, symbolic variables can be decision variables of an optimization problem. The following syntax creates polynomials which have symbolic variables as coefficients. In all of the following syntaxes, the first argument corresponds to the display name (resp., prefix) for the symbolic coefficients. See `Casadi's SX symbolics <https://web.casadi.org/docs/#the-sx-symbolics>`_ for details.

.. code-block:: matlab

    casos.PS.sym('c',w)

creates a scalar polynomial with symbolic coefficients and monomials in ``w``, where ``w`` must be a vector of monomials.

.. code-block:: matlab

    casos.PS.sym('c',w,[m n])
    casos.PS.sym('c',w,n)

creates a ``m × n`` matrix (resp., a square matrix with length ``n``) of polynomials with symbolic coefficients and monomials in ``w``, where ``w`` must be a vector of monomials.

.. code-block:: matlab

    casos.PS.sym('c',[m n])
    casos.PS.sym('c',n)

creates a ``m × n`` matrix (resp., a square matrix with length ``n``) of polynomials of degree zero; essentially, this is a symbolic matrix similar to ``casadi.SX``.

.. code-block:: matlab

    casos.PS.sym(...,'gram')

where ``...`` denotes any of the syntaxes above, creates a scalar or matrix polynomial in Gram form, that is, with entries ``p = z'*Q*z``, where ``z`` is the vector of monomials in ``w`` and ``Q`` is a quadratic symbolic matrix.

.. note::

   The syntax ``casos.PS.sym('c',w,...)`` creates polynomials with monomials *in* ``w`` and is therefore equivalent to ``casos.PS.sym('c',monomials(w),...)``. In consequence, all of the following syntaxes yield the same result:

.. code-block:: matlab
   
   casos.PS.sym('c',[1;w])
   casos.PS.sym('c',[w;1])
   casos.PS.sym('c',[1;w;w])


.. note::
   We say that a polynomial is *symbolic* if and only if all of its (nonzero) coefficients are symbols in the sense of Casadi. Except for the Gram form, all polynomials created with the syntaxes above are symbolic but the result of the notation 

.. code-block:: matlab
   
   casos.PS.sym('c',[1 2])*[x;x]

with ``x = casos.PS('x')`` would only be a symbolic *expression*. The same is true for the Gram form syntax because of the symmetries in the Gram matrix expression. The queries ``is_symbolic``, ``is_symexpr``, and ``is_symgram`` check whether a polynomial is a symbolic polynomial, a symbolic expression, or a symbolic Gram form, respectively.


Class Documentation
-------------------
.. automodule:: casos

.. autoclass:: PS
   :members:
   :exclude-members: casadi.SX, casadi.DM
