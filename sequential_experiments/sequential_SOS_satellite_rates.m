
% clear
close all
clc
profile off

import casos.toolboxes.sosopt.cleanpoly
import casos.toolboxes.sosopt.plinearize
import casos.toolboxes.sosopt.pcontour
import casos.toolboxes.sosopt.pcontour3

% system states
x = casos.Indeterminates('x',3,1);
u = casos.Indeterminates('u',3,1);

% Casini Parameter
J = diag([8970;9230;3830]);

% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

omega_max = 2*pi/180;

% system dynamics
f = -inv(J)*skew(x(1:3))*J*x(1:3);

gx = inv(J);

[A,B] = plinearize(f+gx*u ,x , u);

[K0,P] = lqr(A,B,eye(3)*0.1,eye(3)*0.01);

K = -K0*x;


%% scaled closed-loop dynamics and scaled Lyapunov function
D = diag([1/(omega_max*2),1/(omega_max*2),1/(omega_max*2)]);

% substitute control law and clean up
fc = cleanpoly(subs(f + gx*u,u,K), 1e-10);

% scaling
fc = D*subs(fc,x,D^(-1)*x);


Vinit = (D^(-1)*x)'*P*(D^(-1)*x);

figure(1)
pcontour(subs(Vinit,x(3),0),1,[-omega_max omega_max -omega_max omega_max])
hold on

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));%) monomials(x,2:4));

% SOS multiplier
s2 = casos.PS.sym('s2',monomials(x,2));

s3 = casos.PS.sym('s3',monomials(x,0:2));

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
b = casos.PS.sym('b');

% options
opts = struct('sossol','mosek');


g0 = 0;
n0 = 2;
for k = 1:3

    g0  = g0 + (x(k)^2/omega_max^2)^(n0/2) ;

end
g0 = g0-1;
 

g0 = subs(g0,x,D^(-1)*x);

pcontour(subs(g0,x(3),0),0,[-omega_max omega_max -omega_max omega_max],'k')

cost = dot(g0 - (V), g0 - (V)) ;

%% setup solver
sos1 = struct('x',[V;s2],...
              'f',cost, ...
              'p',[]);

sos1.('g') = [s2; 
              V - l; 
              s2*(V - 1) - nabla(V,x)*fc - l];

% states + constraint are SOS cones


opts.Kx      = struct('lin', length(sos1.x));
opts.Kc      = struct('sos', length(sos1.g));
opts.verbose = 1;
opts.sossol_options.sdpsol_options.error_on_fail = 0;



tic
S1 = casos.nlsossol('S1','sequential',sos1,opts);
toc

sol1 = S1('x0' ,[Vinit; (x'*x)]);

Vsol = subs(sol1.x(1),x,D*x);


figure(1)
pcontour(subs(Vsol,x(3),0),1,[-omega_max omega_max -omega_max omega_max],'r')