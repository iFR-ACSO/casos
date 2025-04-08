%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       region-of-attraction for the longitudinal motion 
%                       of the Nasa Generic Transport Model. To increase
%                       the size of the sublevel set we try to minimize the
%                       squared distance to a defined set. Additionally, we
%                       synthesis a linear control law at the same time.
%                       State and control constraints are considered.
%
%   Reference: Modified problem from:
%              Chakraborty, Abhijit and Seiler, Peter and Balas, Gary J.,
%              Nonlinear region of attraction analysis for flight control 
%              verification and validation, Control Engineering Practice,
%              2011, doi: 10.1016/j.conengprac.2010.12.001
%           
%
%--------------------------------------------------------------------------

import casos.toolboxes.sosopt.cleanpoly

%% system states
x = casos.Indeterminates('x',4,1);
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


%% Trim condition
v0      = 45;      % m/s
alpha0  = 0.04924; % rad
q0      = 0;       % rad/s
theta0  = 0.04924; % radg

eta0    = 0.04892; % rad
delta0  = 14.33;   % percent
            

x0 = [v0; alpha0; q0; theta0];
u0 = [eta0; delta0];

% get trim point 
f = @(x,u) [
                    f1(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
                    f2(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
                    f3(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
                    f4(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
];

[x0, u0] = findtrim(f,x0, u0);

% set up dynamic system (might take some time)
f = f(x+x0,u+u0);

% scale dynamics 
d = ([
          convvel(20, 'm/s', 'm/s')  %range.tas.lebesgue.get('ft/s')
          convang(20, 'deg', 'rad')  %range.gamma.lebesgue.get('rad')
          convang(50, 'deg', 'rad')  %range.qhat.lebesgue.get('rad')
          convang(20, 'deg', 'rad')  %range.alpha.lebesgue.get('rad')
]);



D = diag(d)^-1;

f = D*subs(f, x, D^-1*x);

% remove small terms
f = cleanpoly(f,1e-6,1:5);

% control constraints (shifted by trim)
d2r = pi/180;
uM =   [10*d2r - eta0; 30 - delta0];
um =   [-(10*d2r + eta0); -(10 + delta0) ];

% desired set we want to reach (already scaled)
G0 = [3.95382671347060	-0.230032507976836	 0.0316965275615692	 0.292992065909380
    -0.230032507976836	 0.904764915483640	-0.161191428789579	-1.32594376986430
     0.0316965275615692	-0.161191428789579	 0.0344002648214491	 0.242292666551058
     0.292992065909380	-1.32594376986430	 0.242292666551058	 2.02114797577027];

g = (x'*G0*x) - 1;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2:4));

% SOS multiplier
s2    = casos.PS.sym('s2',monomials(x,4));
s3    = casos.PS.sym('s3',monomials(x,0),[2 1]);
s4    = casos.PS.sym('s4',monomials(x,0),[2 1]);
s5    = casos.PS.sym('s4',monomials(x,0:2));

kappa1 = casos.PS.sym('k',monomials([x(1)],1));
kappa2 = casos.PS.sym('k',monomials([x(3)],1));

kappa = [kappa1;kappa2];

b     = casos.PS.sym('b');

% enforce positivity
l = 1e-6*(x'*x);


%% setup solver
% options
opts = struct('sossol','mosek');
opts.verbose  = 1;
opts.max_iter = 100;

% cost
cost = dot(g-(V-b),g-(V-b)) ;

sos = struct('x',[V; s2;s3;s4;s5;kappa;b],... % decision variables
              'f',cost, ...                   % cost
              'p',[]);                        % parameter

% SOS constraints
sos.('g') = [s2;
             s3;
             s4;
             s5;
             V-l; 
             s2*(V-b)-nabla(V,x)*subs(f,u,kappa)-l;
             s3*(V-b) + kappa - um;
             s4*(V-b) + uM- kappa;
             s5*(V-b) - g
             ];

% states + constraint cones
opts.Kx      = struct('lin', length(sos.x));
opts.Kc      = struct('sos', length(sos.g));

% setup solver
S = casos.nlsossol('S1','filter-linesearch',sos,opts);

%% solve problem

% initial guess
V0 = g;
s20 = (x'*x)^2;
s30 = [1;1];
s40 = [1;1];
s50 = x'*x;
K0  = [1*x(1);1*x(3)];
b0  = 1;

x0 = casos.PD([V0;s20;s30;s40;s50;K0;b0]);

% solve problem
sol = S('x0' ,x0);

casos.postProcessSolver(S,true);

S.stats

S.stats.single_iterations{end}.conic

%% plotting
figure

% re-scale solution
xd = D*x;

Vfun = to_function(subs(sol.x(1),x,xd));
gfun = to_function(subs(g,x,xd));

fcontour(@(x2,x3) full(Vfun(0,x2,x3,0) ), [-1 1 -4 4 ], 'b-', 'LevelList', full(sol.x(end)))
hold on
fcontour(@(x2,x3)  full(gfun(0,x2,x3,0) ), [-1 1 -4 4 ], 'r-', 'LevelList', 0)
hold off
legend('Lyapunov function','Safe set function')
