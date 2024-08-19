%---------------------------------------------------------------------
% sosoptdemo1
%
% Demonstration of SOSOPT function for solving sos optimizations.
% This example uses SOSOPT to solve an SOS feasibility problem for
% a degree four polynomial in two variables.
% See SOSOPT help for more details on the function syntax.
%
% Reference: 
% S. Prajna, A. Papachristodoulou, P. Seiler, and P.A. Parrilo,
% SOSTOOLS: SOSDEMO1, Section 3.1 of SOSTOOLS Version 2.00
% User's Guide, 2004.
%---------------------------------------------------------------------

% Define polynomial variables and polynomial
pvar x1 x2;          
p = 2*x1^4 + 2*x1^3*x2 - x1^2*x2^2 + 5*x2^4;

% Set options for sosopt
opts = sosoptions;
%opts.form = 'image';
opts.form = 'kernel';
%opts.solver = 'csdp';
%opts.solver = 'dsdp';
%opts.solver = 'sdpt3';
%opts.solver = 'sdplr';
%opts.solver = 'sdpam'; %opts.solveropts.epsilonStar = 1e-9;

% Use SOSOPT to check if p is an SOS
pconstr = p>=0;
[info,dopt,sossol] = sosopt(pconstr,opts);

% Check for feasibility of solution
% If p(x1,x2) is SOS then the optimization will return a Gram Matrix 
% decomposition in sossol
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible: p is a sum of squares\n');
    z = sossol.z;
    Q = sossol.Q;
else
    fprintf('\nSOS Optimization is not feasible: p is NOT a sum of squares\n');    
    return
end

% Verify the Gram Marix Decomposition in sossol
e = p - z'*Q*z;
fprintf('\nMax magnitude coefficient of p-z''*Q*z is:')
disp(full(max(abs(e.coefficient))))
fprintf('Minimum eigenvalue of Q is:')
disp(min(eig(Q)));

