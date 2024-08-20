%---------------------------------------------------------------------
% sosoptdemo4
%
% Demonstration of SOSOPT function for solving sos optimizations.
% This example uses SOSOPT to verify the copositivity of a matrix.
% See SOSOPT help for more details on the function syntax.
%
% Reference: 
% S. Prajna, A. Papachristodoulou, P. Seiler, and P.A. Parrilo,
% SOSTOOLS: SOSDEMO4, Section 3.4 of SOSTOOLS Version 2.00
% User's Guide, 2004.
%---------------------------------------------------------------------

% Define polynomial variables
pvar x1 x2 x3 x4 x5;
x = [x1; x2; x3; x4; x5];

% Matrix to be tested
J = [1 -1  1  1 -1;
    -1  1 -1  1  1;
     1 -1  1 -1  1;
     1  1 -1  1 -1;
    -1  1  1 -1  1];


% Define SOS constraint: r(x)*J(x) >=0
v = x.^2;
Jx = v'*J*v;
r = x'*x;
p = r*Jx;
pconstr = p>=0;

% Set options for SOSOPT
% When form='image', the solver runs into numerical
% problems but the solution (which can be obtained by
% turning of the feasibility check) is still feasible.
% Alternatively, the numerical scaling can be turned on.

opts = sosoptions;
%opts.form = 'image'; opts.checkfeas = 'off';
%opts.form = 'image'; opts.scaling = 'on';
%opts.form = 'kernel';
%opts.solver = 'sdpt3';
%opts.solver = 'csdp';
%opts.solver = 'sdplr';

% Call sosopt to test feasibility
[info,dopt,sossol] = sosopt(pconstr,opts);

% Program is feasible. The matrix J is copositive.
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible:\n');
else
    fprintf('\nSOS Optimization is not feasible. \n');    
    return
end

% Verify the Gram Marix Decompositions in sossol
z = sossol.z;
Q = sossol.Q;
e = p - z'*Q*z;
fprintf('\nMax magnitude coefficient of p-z''*Q*z is:')
disp(full(max(abs(e.coefficient))))
fprintf('Minimum eigenvalue of Q is:')
disp(min(eig(Q)));





