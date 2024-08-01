
close all
clc
profile off

import casos.toolboxes.sosopt.cleanpoly
import casos.toolboxes.sosopt.plinearize
import casos.toolboxes.sosopt.pcontour
import casos.toolboxes.sosopt.pcontour3

% system states
x = casos.Indeterminates('x',7,1);
u = casos.Indeterminates('u',3,1);

x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4);
x5 = x(5);
x6 = x(6);
x7 = x(7);

% Casini Parameter
J = diag([8970;9230;3830]);

% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

omega_max = 3*pi/180;


% dynamics
x_dot =  [-inv(J)*skew(x(1:3))*J*x(1:3) + inv(J)*u;                  % omega_dot
           1/2*[-skew(x(1:3)) ,1*x(1:3) ; -1*x(1:3)' ,0]*1*x(4:7)];  % quaternion kinematics

% x0    = [0;0;0;1;0;0;0];
% u0    = zeros(3,1);
% [A,B] = plinearize(x_dot ,x , u,x0,u0);
% 
% [K0,P] = lqr(A,B,diag([0.1,0.1,0.1,1,1,1,1]),eye(3)*0.01);

Kp = 1/2*J;
Kd = diag([1,0.5,2]);

uc = -Kd*x(1:3) - Kp*x(5:7);

f = subs(x_dot,u,uc);

% x0    = [0;0;0;0;0;0;1];
% u0    = zeros(3,1);
% [A,B] = plinearize(f ,x , u,x0,u0);


% lyap(A,eye(7))

%% scaled closed-loop dynamics and scaled Lyapunov function
% D = diag([1/(omega_max*2),1/(omega_max*2),1/(omega_max*2),1,1,1,1]);
D = eye(7);
% substitute control law and clean up
fc = cleanpoly(f, 1e-10);

% scaling
fc = D*subs(fc,x,D^(-1)*x);


% Vinit = (D^(-1)*x)'*P*(D^(-1)*x);

% figure(1)
% pcontour(subs(subs(Vinit,x,D*x),x(3:end),zeros(4,1)),1,[-omega_max omega_max -omega_max omega_max],'b')
% hold on

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s2 = casos.PS.sym('s2',monomials(x,2));
s3 = casos.PS.sym('s3',monomials(x,2));

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

pcontour(subs(subs(g0,x,D*x),x(3:end),zeros(5,1)),0,[-omega_max omega_max -omega_max omega_max]*2,'k')

cost = dot(g0 - (V-1), g0 - (V-1)) ;

%% setup solver
sos = struct('x',[V;s2;s3],...
              'f',cost, ...
              'p',[]);

sos.('g') = [s2; 
             s3;
             V - l; 
             s2*(V - 1) - nabla(V,x)*fc - l;
             s3*(V - 1) + 1 - x(4:end)'*x(4:end)
             ];

% states + constraint are SOS cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));
opts.verbose = 1;

opts.sossol_options.sdpsol_options.error_on_fail = 0;

% profile on
buildTime_in = tic;
    solver_Satellite6D  = casos.nlsossol('S','sequential',sos,opts);
buildtime = toc(buildTime_in);

% profile viewer
% solve problem
sol = solver_Satellite6D ('x0' ,[(x'*x); (x'*x); (x'*x)]);
disp(['Solver buildtime: ' num2str(buildtime), ' s'])

% plot solver statistics
plotSolverStats(solver_Satellite6D.stats);

% plot solution
Vsol = subs(sol.x(1),x,D*x);


figure(1)
pcontour(subs(Vsol,x(3:end),zeros(5,1)),1,[-omega_max omega_max -omega_max omega_max]*2,'r')
