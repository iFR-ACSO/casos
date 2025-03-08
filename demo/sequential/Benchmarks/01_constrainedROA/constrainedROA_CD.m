%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for the longitudinal motion 
%                       of the Nasa Generic Transport Model. To increase
%                       the size of the sublevel set we try to minimize the
%                       squared distance to a defined set. Additionally, we
%                       synthesis a control law at the same time.
%
%   Reference: Modified problem from:
%              Chakraborty, Abhijit and Seiler, Peter and Balas, Gary J.,
%              Nonlinear region of attraction analysis for flight control 
%              verification and validation, Control Engineering Practice,
%              2011, doi: 10.1016/j.conengprac.2010.12.001
%           
%
%--------------------------------------------------------------------------

clc
clear

%import casos.toolboxes.sosopt.plinearize
%import casos.toolboxes.sosopt.pcontour
import casos.toolboxes.sosopt.cleanpoly
%import casos.toolboxes.sosopt.pcontour3

% system states
x = casos.Indeterminates('x',6,1); % we have two additional states now
u = casos.Indeterminates('u',2,1);

%% Polynomial Dynamics
f1 = @(V,alpha,q,theta,eta,F) +1.233e-8*V.^4.*q.^2         + 4.853e-9.*alpha.^3.*F.^3    + 3.705e-5*V.^3.*alpha.*q      ...
                             - 2.184e-6*V.^3.*q.^2         + 2.203e-2*V.^2.*alpha.^3     - 2.836e-6*alpha.^3.*F.^2      ...
                             + 3.885e-7*alpha.^2.*F.^3     - 1.069e-6*V.^3.*q            - 4.517e-2*V.^2.*alpha.^2      ...
                             - 2.140e-3*V.^2.*alpha.*eta   - 3.282e-3*V.^2.*alpha.*q     - 8.901e-4*V.^2.*eta.^2        ...
                             + 9.677e-5*V.^2.*q.^2         - 2.037e-4*alpha.^3.*F        - 2.270e-4*alpha.^2.*F.^2      ...
                             - 2.912e-8*alpha.*F.^3        + 1.591e-3*V.^2.*alpha        - 4.077e-4*V.^2.*eta           ...
                             + 9.475e-5*V.^2.*q            - 1.637   *alpha.^3           - 1.631e-2*alpha.^2.*F         ...
                             + 4.903   *alpha.^2.*theta    - 4.903   *alpha.*theta.^2    + 1.702e-5*alpha.*F.^2         ...
                             - 7.771e-7*F.^3               + 1.634   *theta.^3           - 4.319e-4*V.^2                ...
                             - 2.142e-1*alpha.^2           + 1.222e-3*alpha.*F           + 4.541e-4*F.^2                ...
                             + 9.823   *alpha              + 3.261e-2*F                  - 9.807.*theta      + 4.284e-1;

                      
% change of angle of attack
f2 = @(V,alpha,q,theta,eta,F) -3.709e-11*V.^5.*   q.^2     + 6.869e-11*V.*alpha.^3.*F.^3 + 7.957e-10*V.^4.*alpha.*q     ...
                             + 9.860e-9 *V.^4.*   q.^2     + 1.694e-5 *V.^3.*alpha.^3    - 4.015e-8 *V.*alpha.^3.*F.^2  ...
                             - 7.722e-12*V.*alpha.^2.*F.^3 - 6.086e-9 *alpha.^3.*F.^3    - 2.013e-8 *V.^4.*q            ...
                             - 5.180e-5 *V.^3.*alpha.^2    - 2.720e-6 *V.^3.*alpha.*eta  - 1.410e-7 *V.^3.*alpha.*q     ...
                             + 7.352e-7 *V.^3.*eta.^2      - 8.736e-7 *V.^3.*q.^2        - 1.501e-3 *V.^2.*alpha.^3     ...
                             - 2.883e-6 *V.*alpha.^3.*F    + 4.513e-9 *V.*alpha.^2.*F.^2 - 4.121e-10*V.*alpha.*F.^3     ...
                             + 3.557e-6 *alpha.^3.*F.^2    + 6.841e-10*alpha.^2.*F.^3    + 4.151e-5 *V.^3.*alpha        ...
                             + 3.648e-6 *V.^3.*eta         + 3.566e-6 *V.^3.*q           + 6.246e-6 *V.^2.*alpha.*q     ...
                             + 4.589e-3 *V.^2.*alpha.^2                                  - 6.514e-5 *V.^2.*eta.^2       ...
                             + 2.580e-5 *V.^2.*q.^2        - 3.787e-5 *V.*alpha.^3       + 3.241e-7 *V.*alpha.^2.*F     ...
                             + 2.409e-7 *V.*alpha.*F.^2    + 1.544e-11*V.*F.^3           + 2.554e-4 *alpha.^3.*F        ...
                             - 3.998e-7 *alpha.^2.*F.^2    + 3.651e-8 *alpha.*F.^3       + 4.716e-7 *V.^3               ...
                             - 3.677e-3 *V.^2*alpha        - 3.231e-4 *V.^2.*eta         - 1.579e-4 *V.^2.*q            ...
                             + 2.605e-3 *V.*alpha.^2       + 1.730e-5 *V.*alpha.*F       - 5.201e-3 *V.*alpha.*theta    ...
                             - 9.026e-9 *V.*F.^2           + 2.601e-3 *V.*theta.^2       + 3.355e-3 *alpha.^3           ...
                             - 2.872e-5 *alpha.^2.*F       - 2.134e-5 *alpha.*F.^2       - 1.368e-9 *F.^3               ...
                             - 4.178e-5 *V.^2              + 2.272e-4 *V.*alpha          - 6.483e-7 *V.*F               ...
                             - 2.308e-1 *alpha.^2          - 1.532e-3 *alpha.*F          + 4.608e-1 *alpha.*theta       ...
                             - 2.304e-1 *theta.^2          + 7.997e-7 *F.^2              - 5.210e-3 *V                  ...
                             - 2.013e-2 *alpha             + 5.744e-5 *F          + q    + 4.616e-1;
                         
% change of pitch rate                         
f3 = @(V,alpha,q,theta,eta,F) -6.573e-9*V.^5.*    q.^3     + 1.747e-6*V.^4.*q.^3         - 1.548e-4*V.^3.*q.^3          ...
                             - 3.569e-3*V.^2.*alpha.^3     + 4.571e-3*V.^2.*q.^3         + 4.953e-5*V.^3.*q             ...
                             + 9.596e-3*V.^2.*alpha.^2     + 2.049e-2*V.^2.*alpha.*eta   - 2.431e-2*V.^2.*alpha         ...
                             - 3.063e-2*V.^2.*eta          - 4.388e-3*V.^2.*q            - 2.594e-7*F.^3                ...
                             + 2.461e-3*V.^2               + 1.516e-4*F.^2               + 1.089e-2*F        + 1.430e-1;                         


% change of pitch angle
f4 = @(V,alpha,q,theta,eta,F)  q;

%% Trim condition (original system)
v0      = 45;       % m/s
alpha0  = 0.04924; % rad
q0      = 0;       % rad/s
theta0  = 0.04924; % radg

eta0    = 0.04892; % rad
delta0  = 14.33;   % percent
            

% trim is now augmented by two additional states i.e. trim values
x0 = [v0; alpha0; q0; theta0;eta0; delta0];
u0 = zeros(2,1);

% get trim point 
f = @(x,u) [
                    f1(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:))
                    f2(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:))
                    f3(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:))
                    f4(x(1,:),x(2,:),x(3,:),x(4,:),x(5,:),x(6,:))
                    u
];

[x0, u0] = findtrim(f,x0, u0);

% set up dynamic system
f = f(x+x0,u+u0);

d = ([
          convvel(20, 'm/s', 'm/s')  % range.tas.lebesgue.get('ft/s')
          convang(20, 'deg', 'rad')  % range.gamma.lebesgue.get('rad')
          convang(50, 'deg', 'rad')  % range.qhat.lebesgue.get('rad')
          convang(20, 'deg', 'rad')  % range.alpha.lebesgue.get('rad')
          1                          % no scale for control
          1                          % no scale for control
]);


% use scaled dynamics to compute initial guess for lyapunov function
D = diag(d)^-1;

f = D*subs(f, x, D^-1*x);

f = cleanpoly(f,1e-6,1:5);

% get an initial Lyapunov function for the augmented system
A = full(subs(nabla(f,x),[x;u],[x0;u0]));
B = full(subs(nabla(f,u),[x;u],[x0;u0]));

[K0,P] = lqr(full(A),full(B),eye(6),eye(2));

% initial controller for augmented system
K = -K0*x;

% initial Lyapunov function
Vinit = x'*P*x;


% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));

% SOS multiplier
s1    = casos.PS.sym('s2',monomials(x,4));

% setup linear control laws
k1 = casos.PS.sym('k1',monomials(x(1)),1) ;
k2 = casos.PS.sym('k2',monomials(x(3)),1) ;

kappa = [k1;k2];

s2    = casos.PS.sym('s3',monomials(x,0),[2 1]);
s3    = casos.PS.sym('s4',monomials(x,0),[2 1]);
s4    = casos.PS.sym('s4',monomials(x,0:2));

b    = casos.PS.sym('b');

% enforce positivity
l = 1e-6*(x'*x);

% control constraints (shifted by trim)
d2r = pi/180;
uM =   [  10*d2r - eta0;    30 - delta0];
um =   [-(10*d2r + eta0); -(10 + delta0) ];

% desired set we want to reach (already scaled)
G0 = [3.95382671347060	-0.230032507976836	 0.0316965275615692	 0.292992065909380
    -0.230032507976836	 0.904764915483640	-0.161191428789579	-1.32594376986430
     0.0316965275615692	-0.161191428789579	 0.0344002648214491	 0.242292666551058
     0.292992065909380	-1.32594376986430	 0.242292666551058	 2.02114797577027];

% level set can be used to increase size (but other functions are also
% possible); since have the augmented states, we only consider the original
% states here
g = (x(1:4)'*G0*x(1:4)) - 1;

%% setup solver
bval = 0.01;

% solve for s2 and kappa
opts = struct('sossol','mosek');
opts.verbose = 1;
opts.error_on_fail = 0;
% opts.conf_interval = [0 1];
kappa = K;
sos = struct('x',[s1;s2;s3;s4],...
             'f',-b,...
              'p',V);

sos.('g') = [
             s1; 
             % s2;
             % s3;
             % s4;
             s1*(V-b)-nabla(V,x)*subs(f,u,kappa)-l;
             % s2*(V-b) + kappa - um;
             % s3*(V-b) + uM - kappa;
             % s4*(V-b) - g
             ];

opts.Kx      = struct('lin', 2 + 2*length(u));
opts.Kc      = struct('sos', length(sos.g));

% solver for multiplier and control law step
solver_GTM_syn1 = casos.qcsossol('S1','bisection',sos,opts);

opts = [];
cost = dot(g-(subs(V,x(5:6),zeros(2,1))-b),g-(subs(V,x(5:6),zeros(2,1))-b));

% cost = dot(g-(V-b),g-(V-b));

sos = struct('x',[V;b],...
              'f',cost , ...
              'p',[s1;s2;s3;s4]);

sos.('g') = [V-l; 
             s1*(V-b)-nabla(V,x)*subs(f,u,kappa)-l;
             s2*(V-b) + kappa - um;
             s3*(V-b) + uM - kappa;
             s4*(V-b) - g];

% states + constraint are SOS cones
opts.Kx      = struct('lin', 2);
opts.Kc      = struct('sos', length(sos.g));

% solver for V and beta
solver_GTM_syn2 = casos.sossol('S2','mosek',sos,opts);

Vval = Vinit;

fval_old = [];

for iter = 1:100

% solve for multiplier s2 and control law
sol1 = solver_GTM_syn1('p',Vval);  

if ~strcmp(solver_GTM_syn1.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
    fprintf('Not feasible in multiplier step in iteration %d\n',iter)
    break
end

% solve for V and beta
sol2 = solver_GTM_syn2('p',sol1.x);

if strcmp(solver_GTM_syn2.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS')
   % extract solution
   Vval = sol2.x(1);
    bval = sol2.x(2);

else
    fprintf('Not feasible in Lyapunov step in iteration %d\n',iter)
    break
end

% check convergence of cost function
 if ~isempty(fval_old)
    if abs(full(sol2.f-fval_old)) <= 1e-3
        break
    else
        fval_old = sol2.f;
    end
else
    fval_old = sol2.f;
 end

% show progress 
fprintf('Iteration %d: f = %g, b = %g.\n',iter,full(sol2.f),full(bval));

end


%% plot solver statistics
% plotSolverStats(solver_GTM_syn.stats);

% %% plotting
% import casos.toolboxes.sosopt.pcontour
% 
% xD = D*x;
% V = subs(sol.x(1),[x(1);x(4)],zeros(2,1));
% V = subs(V,[x(2);x(3)],xD(2:3));
% 
% 
% g = subs(g,[x(1);x(4)],zeros(2,1));
% g = subs(g,[x(2);x(3)],xD(2:3));
% 
% figure()
% clf
% pcontour(V, 1, [-4 4 -4 4], 'b-');
