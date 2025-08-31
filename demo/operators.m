% This script demonstrate the usability of linear operators.

x = casos.Indeterminates('x');

% we define a symbolic univariate polynomial
v = casos.PS.sym('c',monomials(x,0:5));

% and take the derivative and antiderivative
% these are symbolic expressions, depending on the symbolic polynomial
dv = nabla(v,x);
iv  =   int(v,x);

% we define derivation and integration operators for v as
D = jacobian(dv,v);
I = jacobian(iv,v);

% these operators can be applied to any other polynomial
% (note: the monomials should not change)
f = casos.PS(sparsity(v),randn(6,1));

disp('The derivative of f is')
disp(evaluate(D,f))
disp(' ')

disp('And the antiderivative of f is')
disp(evaluate(I,f))
disp(' ')

% operators can be added to form a new operator
disp('The sum of the derivative and the antiderivative of is')
disp(evaluate(D + I,f))
disp(' ')

% as well as 'chained' together via composition
% (here, the derivative partially reverses the antiderivative)
DI = compose(D,I);
ID = compose(I,D);

disp('The deriviative of the antiderivative of is')
disp(evaluate(DI,f))
disp(' ')

disp('Whereas the antideriviative of the derivative of is')
disp(evaluate(ID,f))
disp(' ')