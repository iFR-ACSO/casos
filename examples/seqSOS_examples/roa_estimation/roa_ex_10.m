%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for a satellite. We consider a
%                       linear control law. The dynamics consists of the
%                       rotational dynamics and attitude kinematics
%                       (Modified Rodriques Parameter) with respect to the
%                       inertial frame. We also consider constraints on the
%                       angular rates.
%
%   Reference: custom
%           
%
%--------------------------------------------------------------------------

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

[A,B] = plinearize(f + gx*u ,x , u);

[K0,P] = lqr(A,B,eye(3)*0.1,eye(3)*0.01);

K = -K0*x;


%% scaled closed-loop dynamics and scaled Lyapunov function
D = diag([1/(omega_max),1/(omega_max),1/(omega_max)]);

% substitute control law and clean up
fc = cleanpoly(subs(f + gx*u, u, K), 1e-10);

% scaling
fc = D*subs(fc,x,D^(-1)*x);

Vinit = (D^(-1)*x)'*P*(D^(-1)*x);

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s2 = casos.PS.sym('s2',monomials(x,2));
s3 = casos.PS.sym('s3',monomials(x,2));

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
b = casos.PS.sym('b');

% options
opts = struct('sossol','mosek');

% pre-computed safe set
g0 =  673551*x(1)^4 + 673551*x(2)^4 + 673551*x(3)^4 - 0.999933;

% de-scale input
g0 = subs(g0,x,D^(-1)*x);

pcontour(subs(g0,x(3),0),0,[-omega_max omega_max -omega_max omega_max],'k')

cost = dot(g0 - V, g0 - V) ;

%% setup solver
sos = struct('x',[V;s2;s3],...
              'f',cost, ...
              'p',[]);


sos.('g') = [s2;
              s3; 
              V - l; 
              s2*(V - 1) - nabla(V,x)*fc - l;
              s3*(V - 1) - g0];

% states + constraint are SOS cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));

opts.verbose = 1;
opts.indeterminates = casos.PS(x);
opts.sossol_options.sdpsol_options.error_on_fail = 0;


buildTime_in = tic;
    solver_Satellite3D = casos.nlsossol('S','sequential',sos,opts);
buildtime = toc(buildTime_in);

% solve
sol = solver_Satellite3D('x0' ,[Vinit; x'*x;x'*x]);
disp(['Solver buildtime: ' num2str(buildtime), ' s'])

% plot solver statistics
plotSolverStats(solver_Satellite3D.stats);


% Descale/ scale input
Vsol = subs(sol.x(1),x,D*x);
g0   = subs(g0,x,D*x);

figure()
pcontour(subs(Vsol,x(3),0),1,[-omega_max omega_max -omega_max omega_max],'r')
hold on
pcontour(subs(g0,x(3),0),0,[-omega_max omega_max -omega_max omega_max],'k--')