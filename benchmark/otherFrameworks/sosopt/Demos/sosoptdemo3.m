%---------------------------------------------------------------------
% sosoptdemo3
%
% Demonstration of SOSOPT function for solving sos optimizations.
% This example uses SOSOPT to find a lower bound on the global 
% minimum of the Goldstein-Price function.
% See SOSOPT help for more details on the function syntax.
%
% Reference: 
% S. Prajna, A. Papachristodoulou, P. Seiler, and P.A. Parrilo,
% SOSTOOLS: SOSDEMO3, Section 3.3 of SOSTOOLS Version 2.00
% User's Guide, 2004.
%---------------------------------------------------------------------

% Define polynomial variables
pvar x1 x2 gam;
x = [x1;x2];

% Define SOS constraint: f(x) >= gam
% f(x) is the Goldstein-Price function
f1 = x1+x2+1;
f2 = 19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2;
f3 = 2*x1-3*x2;
f4 = 18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2;

f = (1+f1^2*f2)*(30+f3^2*f4);

pconstr = f>=gam;

% Create objective function: 
% maximize gamma <--> minimize (-gamma)
obj = -gam;

% Set options for SOSOPT
opts = sosoptions;
%opts.form = 'image';
%opts.form = 'kernel';
%opts.solver = 'sdpt3';
%opts.solver = 'csdp';
%opts.solver = 'sdplr';
%opts.solver = 'sdpam'; %opts.solveropts.epsilonStar = 1e-9;
%opts.solveropts.print = 'display';
%opts.solveropts.lambdaStar = 1e3;

% Call sosopt to find lower bound on f(x)
[info,dopt,sossol] = sosopt(pconstr,x,obj,opts);

% Get the lower bound
fprintf('\n----------------Results');
if info.feas
    fprintf('\nSOS Optimization is feasible:\n');
    gamsol = subs(gam,dopt)
else
    fprintf('\nSOS Optimization is not feasible. \n');    
    return
end

% Verify the Gram Marix Decomposition in sossol for f(x)-gamsol
z1 = sossol.z;
Q1 = sossol.Q;
e1 = (f-gamsol) - (z1'*Q1*z1);
fprintf('\nMax magnitude coefficient of s(1)-z1''*Q1*z1 is:')
disp(full(max(abs(e1.coefficient))))
fprintf('Minimum eigenvalue of Q1 is:')
disp(min(eig(Q1)));




