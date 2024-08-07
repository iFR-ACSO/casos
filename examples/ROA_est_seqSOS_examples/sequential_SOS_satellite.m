
close all
clc
profile off

import casos.toolboxes.sosopt.cleanpoly
import casos.toolboxes.sosopt.plinearize
import casos.toolboxes.sosopt.pcontour
import casos.toolboxes.sosopt.pcontour3

% system states
x = casos.Indeterminates('x',6,1);
u = casos.Indeterminates('u',3,1);

x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4);
x5 = x(5);
x6 = x(6);

% Casini Parameter
J = diag([8970;9230;3830]);

% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

omega_max = 3*pi/180;

% system dynamics
B = @(sigma) (1-sigma'*sigma)*eye(3)+skew(sigma)+ 2*sigma*sigma';

% dynamics
x_dot =  [-inv(J)*skew(x(1:3))*J*x(1:3) + inv(J)*u; % omega_dot
           1/4*B(x(4:6))*x(1:3)];                   % MRP kinematics


[A,B] = plinearize(x_dot ,x , u);

[K0,P] = lqr(A,B,diag([0.1,0.1,0.1,1,1,1]),eye(3)*0.01);

K = -K0*x;


%% scaled closed-loop dynamics and scaled Lyapunov function
D = diag([1/(omega_max*2),1/(omega_max*2),1/(omega_max*2),1,1,1]);

% substitute control law and clean up
fc = cleanpoly(subs(x_dot,u,K), 1e-10);

% scaling
fc = D*subs(fc,x,D^(-1)*x);


Vinit = (D^(-1)*x)'*P*(D^(-1)*x);

figure(1)
pcontour(subs(subs(Vinit,x,D*x),x(3:end),zeros(4,1)),1,[-omega_max omega_max -omega_max omega_max],'b')
hold on

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));

% SOS multiplier
s2 = casos.PS.sym('s2',monomials(x,2:4));

% enforce positivity
l = 1e-6*(x'*x);

% options
opts = struct('sossol','mosek');


g0 = 0;
n0 = 2;
for k = 1:3

    g0  = g0 + (x(k)^2/omega_max^2)^(n0/2) ;

end
g0 = g0-1;
 

g0 = subs(g0,x,D^(-1)*x);

pcontour(subs(subs(g0,x,D*x),x(3:end),zeros(4,1)),0,[-omega_max omega_max -omega_max omega_max]*2,'k')

cost = dot(g0 - (V-1), g0 - (V-1)) ;

%% setup solver
sos = struct('x',[V;s2],...
              'f',cost, ...
              'p',[]);

sos.('g') = [s2; 
              V - l; 
              s2*(V - 1) - nabla(V,x)*fc - l];

% states + constraint are SOS cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.verbose = 1;

opts.sossol_options.sdpsol_options.error_on_fail = 0;

profile on
buildTime_in = tic;
    solver_Satellite6D  = casos.nlsossol('S','sequential',sos,opts);
buildtime = toc(buildTime_in);

profile viewer
% solve problem
sol = solver_Satellite6D ('x0' ,[Vinit; (x'*x)]);
disp(['Solver buildtime: ' num2str(buildtime), ' s'])

% plot solver statistics
plotSolverStats(solver_Satellite6D.stats);

% plot solution
Vsol = subs(sol.x(1),x,D*x);


figure(1)
pcontour(subs(Vsol,x(3:end),zeros(4,1)),1,[-omega_max omega_max -omega_max omega_max]*2,'r')
