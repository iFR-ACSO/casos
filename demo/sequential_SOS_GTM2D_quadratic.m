% Estimate region of the Generic Transport Model
% See Chakraborty et al. 2011 (CEP) for details.

import casos.toolboxes.sosopt.*

clear
close all
clc


% system states
x = casos.PS('x',2,1);

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

Kq = 0.0698;

%% Trim condition
v0      = 45;       % m/s
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

% find more accurate trim values  
[x0, u0] = findtrim(f,x0, u0);

% short period
f = @(x,u) [
                    f2(x0(1),x(1),x(2),x0(4),u0(1),u0(2))
                    f3(x0(1),x(1),x(2),x0(4),u0(1),u0(2))
];


% set up dynamic system
f = f(x+[x0(2);x0(3)],u0(1));

% scaling
d = ([
          convvel(20, 'm/s', 'm/s')  %range.tas.lebesgue.get('ft/s')
          convang(20, 'deg', 'rad')  %range.gamma.lebesgue.get('rad')
          convang(50, 'deg', 'rad')  %range.qhat.lebesgue.get('rad')
          convang(20, 'deg', 'rad')  %range.alpha.lebesgue.get('rad')
]);


D = diag(d(2:3))^-1;
f = D*subs(f, x, D^-1*x);


f = cleanpoly(f,1e-6,0:5);

% use scaled dynamics to compute initial guess for lyapunov function
A = nabla(f,x);

A0 = double(subs(A,x,zeros(2,1))); 

P = lyap(A0',eye(2));

% initial Lyapunov function
Vinit = x'*P*x;

% polynomial shape
p = x'*x;

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,0));
s2 = casos.PS.sym('s2',monomials(x,2));

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
b = casos.PS.sym('b');

g = Vinit-2; 

cost = dot(g - (V-b),g - (V-b)) - b^2;


%% setup solver
profile on
% options
opts = struct('sossol','mosek');
opts.verbose = 1;

sos1 = struct('x',[V; s2; b],...
              'f',cost, ...
              'p',[]);

sos1.('g') = [s2; 
              V-l; 
              s2*(V-b)-nabla(V,x)*f-l];

% states + constraint are SOS cones
opts.Kx = struct('lin', 3);
opts.Kc = struct('sos', 3);

    
Vlb  = casos.PS(basis(V),-inf);
Vub  = casos.PS(basis(V),+inf);
s1lb = casos.PS(basis(s1),-inf);
s1ub = casos.PS(basis(s1),+inf);
s2lb = casos.PS(basis(s2),-inf);
s2ub = casos.PS(basis(s2),+inf);
glb  = casos.PS(basis(b),-inf);
gub  = casos.PS(basis(b),+inf);

opts.Sequential_Algorithm = 'SQP';

tic
S1 = casos.nlsossol('S1','sequential',sos1,opts);
toc

tic

%% solve
sol1 = S1('x0',[Vinit ;  (x'*x) ; 1], ...
          'lbx',[Vlb;s2lb;glb], ...
          'ubx',[Vub;s2ub;gub]);
toc

profile viewer

% show solution 
% sol1.f
% sol1.x(end)

figure()
pcontour(g,0,[-2 2 -2 2]*3,'k--')
hold on
pcontour(sol1.x(1),double(sol1.x(end)),[-2 2 -2 2]*3)

for k =  1:length(sol1.g)
    solChecker(sol1.g(k));
end

