% SOSDEMO9s --- Matrix SOS decomposition
% Section 4.9 of SOSTOOLS User's Manual

clear; echo on;
syms x1 x2 x3;
vartable = [x1, x2, x3];

% =============================================
% First, initialize the sum of squares program
prog = sosprogram(vartable);   % No decision variables.

% =============================================
% Consider the following candidate sum of squares matrix P(x)
P = [x1^4+x1^2*x2^2+x1^2*x3^2 x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2); 
     x1*x2*x3^2-x1^3*x2-x1*x2*(x2^2+2*x3^2) x1^2*x2^2+x2^2*x3^2+(x2^2+2*x3^2)^2];  

% Test if P(x1,x2,x3) is an SOS matrix and return H so that P = H.'*H
solver_opt.solver = 'sedumi';
[Q,Z,H] = findsos(P,[],solver_opt);

% =============================================
% If program is feasible, P(x1,x2,x3) is an SOS matrix.
echo off;