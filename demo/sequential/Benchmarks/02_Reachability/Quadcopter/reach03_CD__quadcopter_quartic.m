%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       forward reachable set for the planar motion of a 
%                       quadcopter. To increase the size of the sublevel
%                       set we try to minimize the squared distance to 
%                       safe set. Scaled model,
%                       constraints and time horizon from [1] and [2].
%
%   References: 
%          [1] Yin, H., Arcak, M. and Seiler, P. ,
%              Backward Reachability for Polynomial Systems on a Finite Horizon, 
%              IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 66, NO. 12, DECEMBER 2021,
%              doi: 10.1109/TAC.2021.3056611
%          [2] https://github.com/heyinUCB/Backward-Reachability-Analysis-and-Control-Synthesis 
%
%--------------------------------------------------------------------------

clear
clc

% system states 
x = casos.Indeterminates('x',6);

% time
t = casos.Indeterminates('t');

x1 = x(1); % pos. x 
x2 = x(2); % pos. y
x3 = x(3); % vel. x
x4 = x(4); % vel. y
x5 = x(5); % q
x6 = x(6); % theta

% initial guess copter Lyap; see Yin implementation
P = initGuess_copter;

% Polynomial Dynamics
d2r = pi/180;

% 6-state system
g = 9.81;
K = 0.89/1.4;
d0 = 70;
d1 = 17;
n0 = 55;

f = [x3; x4; 0; -g; x6; -d0*x5-d1*x6];

g1 = [0; 0; K*(-0.166*x5^3+x5); K*(-0.498*x5^2+1); 0; 0];
g2 = [0; 0; 0; 0; 0; n0];

% constraints
% rt1 = (1.7 - 1*x1)*(1*x1 + 1.7);
% rt2 = (0.85 - 1*x2)*(1*x2 + 0.85);
% rt3 = (0.8 - 1*x3)*(1*x3 + 0.8);
% rt4 = (1 - 1*x4)*(1*x4 + 1);
% rt5 = (pi/12 - 1*x5)*(1*x5 + pi/12);
% rt6 = (pi/2 - 1*x6)*(1*x6 + pi/2);

% superquadric to cover all constraints
n = 4;
rt = (x(1)^2/1.7^2)^(n/2) + (x(2)^2/0.85^2)^(n/2) + (x(3)^2/0.8^2)^(n/2) + ...
     (x(4)^2/1^2)^(n/2) + (x(5)^2/(pi/12)^2)^(n/2)   + (x(6)^2/(pi/2)^2)^(n/2) - 1;


% initial storage function
G0 = P;
g0 = x'*G0*x - 1;
r0   = rt;

uM = [1.5+g/K;   pi/12];
um = [-1.5+g/K; -pi/12];


% Lyapunov function candidate
V    = casos.PS.sym('v',monomials([x;t],0:4));
K    = casos.PS.sym('k',monomials([x;t],0:2),[2,1]);

% SOS multiplier
s2 = casos.PS.sym('s2',monomials([x;t],0:4));
s3 = casos.PS.sym('s3',monomials([x;t],0:4));
s4 = casos.PS.sym('s4',monomials([x;t],0:2),[length(rt),1]);
s5 = casos.PS.sym('s5',monomials([x;t],0:2),[2,1]);
s6 = casos.PS.sym('s6',monomials([x;t],0:2),[2,1]);
s7 = casos.PS.sym('s7',monomials([x;t],0:2),[2,1]);
s8 = casos.PS.sym('s8',monomials([x;t],0:2),[2,1]);
s9 = casos.PS.sym('s9',monomials([x;t],0:2),[length(rt),1]);
s10 = casos.PS.sym('s10',monomials(x,0:2));
gamma = casos.PS.sym('g');
T = 2;
h = casos.PS(1*t*(T-t*1));

% fixed level set
gamma = 0.1;
% gamma = 0;
%% setup solver
sos1 = struct('x',[K;s2;s3;s4;s5;s6;s7;s8;s9;s10], ... % dec.var
              'p',V);                                   % parameter


% constraint
sos1.('g') = [
                s2;
                s3;
                s4;
                s5;
                s6;
                s7;
                s8;...
                s9;...
                s10;
                s2*rt - s3*h - (nabla(V, t) + nabla(V, x)*( f + g1*K(1) + g2*K(2) ) ) ; ... % partially convexified see Yin paper
                s4.*(V-gamma) - rt      - s9.*h ;...
                s5*(V-gamma)   + uM - K  - s6*h; ...    % convexified see Yin paper
                s7*(V-gamma)  + K - um  - s8*h;         % convexified see Yin paper
                s10*(subs(V,t,0)-gamma) - g0
             ];


% states + constraint are SOS cones
opts.Kx = struct('lin',length(sos1.x));
opts.Kc = struct('sos', length(sos1.g));

% build first solver
measBuildTime1 = tic;
S1 = casos.sossol('S1','mosek',sos1,opts);
buildTime1 = toc(measBuildTime1);


%% setup solver
sos2 = struct('x',V, ...                                  % dec.var
              'f',dot(rt-(V-gamma),rt-(V-gamma)), ...     % cost function for bisection
              'p',[K;s2;s3;s4;s5;s6;s7;s8;s9;s10]); % parameter


% constraint
sos2.('g') = [
                s2*rt - s3*h - (nabla(V, t) + nabla(V, x)*( f + g1*K(1) + g2*K(2) ) ) ; ... % partially convexified see Yin paper
                s4.*(V-gamma) - rt      - s9.*h ;...
                s5*(V-gamma)  + uM - K  - s6*h; ...    % convexified see Yin paper
                s7*(V-gamma) + K - um  - s8*h;         % convexified see Yin paper
                s10*(subs(V,t,0)-gamma) - g0
             ];

opts = [];
% states + constraint are SOS cones
opts.Kx = struct('lin',length(sos2.x));
opts.Kc = struct('sos', length(sos2.g));

% build first solver
measBuildTime2 = tic;
S2 = casos.sossol('S2','mosek',sos2,opts);
buildTime2 = toc(measBuildTime2);

% call solver


import casos.toolboxes.sosopt.*
Vval = x'*P*x;
figure(1)
clf
pcontour(subs(Vval,x(3:6),zeros(4,1)),0,[-1 1 -1 1]*2,'b')
hold on 
pcontour(subs(g0,x(3:6),zeros(4,1)),0,[-1 1 -1 1]*2,'k--')
pcontour(subs(rt,x(3:6),zeros(4,1)),0,[-1 1 -1 1]*2,'k--')


fval_old = [];
startSolve = tic;


load init_solution_seq_quartic.mat
x0 = casos.PD(monoV2,coeffV2);
Vval = x0(1);
for iter = 1:100

% solve for multiplier and control law
sol1 = S1('p',Vval);  

% solve for Lyapunov and minimize cost
sol2 = S2('p',sol1.x);

Vval = sol2.x(1);

% show progress 
fprintf('Iteration %d: f = %g, b = %g.\n',iter,full(sol2.f),full(sol1.f));

% check convergence of cost function
 if ~isempty(fval_old)
    if abs(full(sol2.f-fval_old)) <= 1e-4
        
        break
    else
        fval_old = sol2.f;
    end
else
    fval_old = sol2.f;
 end
end

buildtime = buildTime1 + buildTime2
totalSolveTime = toc(startSolve)
timePerIteration = totalSolveTime/iter

n_con1 = S1.stats.conic.size_A.size(1)
n_dec1 = S1.stats.conic.size_A.size(2)

n_con2 = S2.stats.conic.size_A.size(1)
n_dec2 = S2.stats.conic.size_A.size(2)

%% plotting of sublevel sets
import casos.toolboxes.sosopt.*

Vsol = sol2.x(1)-gamma0;

% initial set
Vval0 = subs(Vsol,t,0);
% reachable set
VvalT = subs(Vsol,t,T);

figure(2)
clf
pcontour(subs(Vval0,x(3:6),zeros(4,1)),0,[-1 1 -1 1]*2,'b')
hold on 
pcontour(subs(VvalT,x(3:6),zeros(4,1)),0,[-1 1 -1 1]*2,'r')
pcontour(subs(g0,x(3:6),zeros(4,1)),0,[-1 1 -1 1]*2,'k')
pcontour(subs(rt,x(3:6),zeros(4,1)),0,[-1 1 -1 1]*2,'k--')
legend('V(0,x)','V(T,x)','X_0','X')
