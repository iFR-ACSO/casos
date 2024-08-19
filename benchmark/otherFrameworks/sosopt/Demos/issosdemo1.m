%---------------------------------------------------------------------
% issosdemo1
%
% Demonstration of the ISSOS function for testing if a polynomial p 
% is a sum of squares.  This example uses ISSOS to construct an
% SOS decomposition for a degree four polynomial in two variables.
% See ISSOS help for more details on the function syntax.
%
% Reference: 
% S. Prajna, A. Papachristodoulou, P. Seiler, and P.A. Parrilo,
% SOSTOOLS: SOSDEMO1, Section 3.1 of SOSTOOLS Version 2.00
% User's Guide, 2004.
%---------------------------------------------------------------------

% Define polynomial variables and polynomial
pvar x1 x2;          
p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4;

% Use SOSOPT to check if p is an SOS
[feas,z,Q,f] = issos(p,opts);

% Check for feasibility of solution
% If p(x1,x2) is SOS then issos will return an SOS decomposition in f
% and a positive semidefinite Gram matrix Q.
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible: p is a sum of squares\n');
    z = sossol.z;
    Q = sossol.Q;
else
    fprintf('\nSOS Optimization is not feasible: p is NOT a sum of squares\n');    
    return
end

% Verify the SOS decomposition is valid
e1 = p - f'*f;
fprintf('\nMax magnitude coefficient of p-f''*f is:')
disp(full(max(abs(e1.coefficient))))

% Verify the Gram Marix Decomposition is valid
e2 = p - z'*Q*z;
fprintf('\nMax magnitude coefficient of p-z''*Q*z is:')
disp(full(max(abs(e2.coefficient))))
fprintf('Minimum eigenvalue of Q is:')
disp(min(eig(Q)));

