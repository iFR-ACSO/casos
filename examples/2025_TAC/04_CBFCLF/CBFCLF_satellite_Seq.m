%% ------------------------------------------------------------------------
%   
%   Supplementary Material for "Infinitesimal-horizon model predictive 
%   control as control barrier and Lyapunov function approach" by 
%   Jan Olucak, Arthur Castello B. de Oliveira, and Torbjørn Cunis
%
%   Short Description: Script to syntheszise the terminal ingredients. The
%                      notation does not fit
%
%
%   Needed software: - CasADi 3.6 
%                    - A C/C++ compiler for .mex file generation (We use a
%                      MS 2022 C/C++ compiler) 
%
%
% ------------------------------------------------------------------------

clear
close all
clc

% system states
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);

%% Hubble telescope parameter
J = diag([3104;7721;7875]);

% simple bounds on rates;
omegaMax1 = 0.5*pi/180;
omegaMax2 = 0.2*pi/180;
omegaMax3 = 0.2*pi/180;

x_low =  [-omegaMax1 -omegaMax2 -omegaMax3]';
x_up  =  [ omegaMax1  omegaMax2  omegaMax3]';


% control constraint; assumption is that the box is inside the full
% torque volume. This is roughly estimated visually.
umin = [-1 -1 -1]'*1.2;
umax = [ 1  1  1]'*1.2;

n = 4;
gu = (u(1)^2/umax(1)^2)^(n/2) + (u(2)^2/umax(2)^2)^(n/2) + (u(3)^2/umax(3)^2)^(n/2) - 1;

Dx   = diag([1/(x_up(1)-x_low(1)),1/(x_up(2)-x_low(2)),1/(x_up(3)-x_low(3)),0.5,0.5,0.5]);

Dxin = inv(Dx);

%% dynamics
% skew-symmetric matrix
cpm = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

% dynamics
B = @(sigma) (1-sigma'*sigma)*eye(3)+ 2*cpm(sigma)+ 2*sigma*sigma';

f =  [-J\cpm(x(1:3))*J*x(1:3) + J\u;
      1/4*B(x(4:6))*x(1:3)]; 

% generate an initial guess for CLF
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

A = full(casos.PD(subs(nabla(f,x),[x;u],[x0;u0])));
B = full(casos.PD(subs(nabla(f,u),[x;u],[x0;u0])));

Q = diag([1, 1, 1, 1,1 ,1]);
R = eye(3)*1;
[~,P0] = lqr(full(A),full(B),Q,R);

% scaled initial guess for terminal penalty (Lyapunov linear system)
Wval = (inv(Dx)*x)'*P0*(inv(Dx)*x);

% scale dynamics
f = Dx*subs(f,[x;u],[Dx\x;u]);

% allowable set (inner-approximation of box constraints via superquadric)
n = 4;
g0 = (x(1)^2/omegaMax1^2)^(n/2) + (x(2)^2/omegaMax2^2)^(n/2) + (x(3)^2/omegaMax3^2)^(n/2) + ...
     (x(4)^2/0.57^2)^(n/2) + (x(5)^2/0.57^2)^(n/2) + (x(6)^2/0.57^2)^(n/2) - 1;

% re-scale input of state constraints
g = subs(g0,x,Dx\x); 


%% setup SOS problem
% terminal set (invariant set)
W  = casos.PS.sym('w',monomials(x,0:4));

% terminal penalty
V  = casos.PS.sym('v',monomials(x,2));

% control law(s)
K1  = casos.PS.sym('k',monomials(x,0),[3 1]);
for j = 1:3
 K1(j) = casos.PS.sym('k',monomials([x(j)]));
end

K2  = casos.PS.sym('k',monomials(x,0),[3 1]);
for j = 1:3
 K2(j) = casos.PS.sym('k',monomials([x(j+3)]));
end

K = K1+K2;

% SOS mulitplier
s1 = casos.PS.sym('s1',monomials(x,0:2));
s2 = casos.PS.sym('s2',monomials(x,2));
s3 = casos.PS.sym('s3',monomials(x,0),[3 1]);
s4 = casos.PS.sym('s4',monomials(x,0),[3 1]);
% s5 = casos.PS.sym('s5',monomials(x,1:2));
s6 = casos.PS.sym('s6',monomials(x,2));

s5 = casos.PS.sym('s5',monomials([x(1)^2 x(2)^2 x(3)^2 x(4)^2 x(5)^2 x(6)^2 x(1)*x(4) x(2)*x(5) x(3)*x(6)]));


% fixed level set of terminal set
b = 0.9;

% options for sequential sos
opts = struct('sossol','mosek');

opts.verbose       = 1;
opts.max_iter      = 100;

cost = dot(g-W,g-W) ;

sos = struct('x', [W; V;K1;K2;s1;s2;s3; s4;s5],... % decision variables
              'f', cost ,...                       % cost function
              'p',[]);                             % parameter

% constraints
sos.('g') = [
             s1
             s3;
             s4;
             s5;
             s1*(W-b)  - g;                         % State constraints
             s2*(W-b)  -  nabla(W,x)*subs(f,u,K);   % set invariance
             s3*g   + K-umin;                   % control constraints
             s4*g   + umax-K;
             s5*g -  nabla(V,x)*subs(f,u,K);  % CLF
             V - 1e-6*(x'*x) % ensure positivity of CLF
             ];

% CLF must be valid on invariant set. Hence, it is also bounded

% states + constraint are linear/SOS cones
opts.Kx = struct('lin', length(sos.x));
opts.Kc = struct('sos', length(sos.g));

% solver setup
S  = casos.nlsossol('S','sequential',sos,opts);

% initial guess for sequential
x0 = casos.PD([ g;  ...
                 Wval;
                -eye(3)*x(1:3); ...         
                -eye(3)*x(4:6); ...
                ones(3,1)*(x'*x);
                ones(3,1)*(x'*x);
                x'*x;x'*x;x'*x]);


% % solve
sol = S('x0',x0);


S.stats
S.stats.single_iterations{end}.SDP_data.size_A

bsol = b;

% re-scale invariant set, terminal penalty and local control law
Wsol_re = subs(sol.x(1),x,Dx*x) - full(casos.PD(bsol));
Vsol_re = subs(sol.x(2),x,Dx*x);
Ksol_re = subs(sol.x(3:5),x,Dx*x) + subs(sol.x(6:8),x,Dx*x);



%% plotting
import casos.toolboxes.sosopt.*

% plot in grad instead of rad; for pcontour the input is given in deg so we  scale the input

% slice for rates
figure(1)
deg2rad = diag([pi/180,pi/180,pi/180 1 1 1]);
clf
pcontour(subs(subs(Wsol_re,x(3:end),zeros(4,1)),x,deg2rad*x),0,[-omegaMax1 omegaMax1 -omegaMax1 omegaMax1]*180/pi,'g')
hold on 
pcontour(subs(subs(g0,x(3:end),zeros(4,1)),x,deg2rad*x),0,[-omegaMax1 omegaMax1 -omegaMax1 omegaMax1]*180/pi,'k--')
legend('Terminal Set','Safe Set')

% 3D slice for Modified rodrigues parameter
figure(2)
deg2rad = diag([pi/180,pi/180,pi/180 1 1 1]);
clf
pcontour3(subs(Wsol_re,x(1:3),zeros(3,1)),0,[-4 4 -4 4 -4 4],'g')
hold on 
legend('Terminal Set')

