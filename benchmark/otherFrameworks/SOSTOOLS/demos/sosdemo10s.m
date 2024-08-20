% SOSDEMO10s --- Set containment
% Section 4.10 of SOSTOOLS User's Manual

clear; echo on;
syms x1 x2;
vartable = [x1 x2];

eps = 1e-6;

% =============================================
% This is the problem data
p = x1^2+x2^2;
gamma = 1;
g0 = [2 0]*[x1;x2];
theta = 1;

% =============================================
% Initialize the sum of squares program
prog = sosprogram(vartable);

% =============================================
% The multiplier
Zmon = monomials(vartable,0:4);
[prog,s] = sospolymatrixvar(prog,Zmon,[1 1]);

% =============================================
% Term to be added to g0
Zmon = monomials(vartable,2:3);
[prog,g1] = sospolymatrixvar(prog,Zmon,[1 1]);

% =============================================
% The expression to satisfy the set containment
Sc = [theta^2-s*(gamma-p) g0+g1; g0+g1 1];

prog = sosmatrixineq(prog,Sc-eps*eye(2));

solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt);

s = sosgetsol(prog,s);
g1 = sosgetsol(prog,g1);

% =============================================
% If program is feasible, { x |((g0+g1) + theta)(theta - (g0+g1)) >=0 } contains { x | p <= gamma }
echo off;