%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction (ROA) for a satellite. 
%                       We consider a linear control law. The dynamics 
%                       consists of the rotational dynamics with respect 
%                       to the inertial frame. Additionally, we synthesis a
%                       control law to increase the size of the ROA. 
%
%   Reference: custom
%           
%
%--------------------------------------------------------------------------


close all
clc
profile off


import casos.toolboxes.sosopt.plinearize
import casos.toolboxes.sosopt.cleanpoly
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
fc = cleanpoly(f + gx*u, 1e-10);

% scaling
fc = D*subs(fc,x,D^(-1)*x);


Vinit = (D^(-1)*x)'*P*(D^(-1)*x);

figure(1)
pcontour(subs(Vinit,x(3),0),1,[-omega_max omega_max -omega_max omega_max])
hold on

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));%) monomials(x,2:4));

% SOS multiplier
s2    = casos.PS.sym('s2',monomials(x,2));
kappa = casos.PS.sym('k',monomials(x,1),[3,1]);

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

pcontour(subs(g0,x(3),0),0,[-omega_max omega_max -omega_max omega_max],'k')

cost = dot(g0 - (V-1), g0 - (V-1)) ;

%% setup solver
sos = struct('x',[V;s2;kappa],...
              'f',cost, ...
              'p',[]);

sos.('g') = [s2; 
              V - l; 
              s2*(V - 1) - nabla(V,x)*subs(fc,u,kappa) - l];

% states + constraint are SOS cones


opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.verbose = 1;
opts.sossol_options.sdpsol_options.error_on_fail = 0;

% build solver
buildTime_in = tic;
    solver_Satellite3D_syn = casos.nlsossol('S','sequential',sos,opts);
buildtime = toc(buildTime_in);

% solve
sol = solver_Satellite3D_syn('x0' ,[x'*x; (x'*x);ones(3,1)]);
disp(['Solver buildtime: ' num2str(buildtime), ' s'])

% plot solver statistics
plotSolverStats(solver_Satellite3D_syn.stats);


% plot solution
Vsol = subs(sol.x(1),x,D*x);
kappa_sol  = subs(sol.x(3:end),x,D*x);

figure(1)
pcontour(subs(Vsol,x(3),0),1,[-omega_max omega_max -omega_max omega_max],'r')
