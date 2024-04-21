% ------------------------------------------------------------------------
%
% Short Description: 
%  Trial to implement contractive sequence of controllable sets.
%
% Date: 04/20/2024
%
% Reference: 
% none
%
% ------------------------------------------------------------------------

import casos.toolboxes.sosopt.plinearize
import casos.toolboxes.sosopt.pcontour

clc
clear
close all

% disp(['--------------------------------------------------------------------' ...
%     '---------'])
%% Define states and load necessary data
% system states
x = casos.PS('x',2,1);
u = casos.PS('u');
t = casos.PS('t');


% system dynamics
f = [x(2);
        (1-x(1)^2)*x(2)-x(1)];

gx = [0;1];

[A,B] = plinearize(f+gx*u,x,u);

[~,P] = lqr(A,B,eye(2),2.5);


umin = -1;
umax =  1;

Omega = x'*P*x - 1; 

% initial Lyapunov functionn
Vval = Omega;


Xs = 3*x(2)^2 + x(1)^2 -1;

%% Define all polynmoials
T = 0.1;
N = 10;


t0 = 0;
T  = T/N;

hT = (t-t0)*(T-t);

% Lyapunov function candidate
V    = casos.PS.sym('v',monomials([x;t],0:3));
Vsym = casos.PS.sym('vsym',monomials([x;t],0:3));
K    = casos.PS.sym('K',monomials(x,0:2));


% SOS multiplier
s1 = casos.PS.sym('s1',monomials([x;t],0:2),'gram');
s2 = casos.PS.sym('s2',monomials([x;t],0:2),'gram');

s3 = casos.PS.sym('s3',monomials([x;t],0:2),'gram');
s4 = casos.PS.sym('s4',monomials([x;t],0:2),'gram');

s51 = casos.PS.sym('s51',monomials([x;t],0:2),'gram');
s52 = casos.PS.sym('s52',monomials([x;t],0:2),'gram');

s61 = casos.PS.sym('s61',monomials([x;t],0:2),'gram');
s62 = casos.PS.sym('s62',monomials([x;t],0:2),'gram');

s7 = casos.PS.sym('s7',monomials(x,0:2),'gram');
s8 = casos.PS.sym('s8',monomials(x,0:2),'gram');
s9 = casos.PS.sym('s9',monomials(x,0:2),'gram');

%% setup solvers for  iter = 0
disp(['--------------------------------------------------------------------' ...
    '---------'])
disp('Setup solver ...')
tic


prob_s_step = struct('x',[K;s1;s2;s3;s4;s51;s52;s61;s62;s7;s8], ...   % decision variables
                     'p',Vsym);                                 % parameter

prob_s_step.('g') = [s1*Vsym - s2*hT - nabla(Vsym,t) - nabla(Vsym,x)*(f + gx*K); 
                     s3*Vsym - s4*hT - Xs;
                    K - umin + s51*Vsym - s61*hT  ; 
                    umax - K + s52*Vsym - s62*hT  ; 
                    s7*subs(Vsym,t,T) - Omega;
                    s8*Omega - subs(Vsym,t,t0)]; % alternative  s8*subs(Vsym,t,T) - subs(Vsym,t,t0)

% states + constraint are SOS cones
opts.Kx = struct('sos', length(prob_s_step.x)-1,'lin',1);
opts.Kc = struct('sos', length(prob_s_step.g));

solver_s_step = casos.sossol('S', ...
                            'mosek', ...
                            prob_s_step, ...
                            opts);



% disp(['... s step succesful after ',  num2str(toc) ' s'])

Ksym  = casos.PS.sym('K',basis(K));
s1sym = casos.PS.sym('s1sym',basis(s1));
s2sym = casos.PS.sym('s2sym',basis(s2));

s3sym = casos.PS.sym('s3sym',basis(s3));
s4sym = casos.PS.sym('s4sym',basis(s4));

s51sym = casos.PS.sym('s51sym',basis(s51));
s52sym = casos.PS.sym('s52sym',basis(s52));

s61sym = casos.PS.sym('s61sym',basis(s61));
s62sym = casos.PS.sym('s62sym',basis(s62));
s7sym  = casos.PS.sym('s7sym',basis(s7));
s8sym  = casos.PS.sym('s8sym',basis(s8));

Vfun0 = to_function(subs(V,t,0));
VfunT = to_function(subs(V,t,T));
gfun  = to_function(Xs);

Nsample = 100;
a = -1;
b = -1;
samples = randn(2,Nsample);
gval = full(casadi.DM(gfun(samples)));

idx = find(gval <= 0);

cost    = 1/Nsample*sum( ( Vfun0(casadi.SX(samples(:,idx)))  ) );


prob_v_step = struct('x',V, ...                % decision variables
                     'f',cost,...
                     'p',[Ksym;s1sym;s2sym;s3sym;s4sym;s51sym;s52sym; ... % parameter
                          s61sym;s62sym;s7sym;s8sym]); 

prob_v_step.('g') = [s1sym*V - s2sym*hT - nabla(V,t) - nabla(V,x)*(f + gx*Ksym); 
                     s3sym*V - s4sym*hT - Xs;
                    Ksym - umin + s51sym*V - s61sym*hT  ; 
                    umax- Ksym  + s52sym*V - s62sym*hT  ; 
                    s7sym*subs(V,t,T) - Omega;
                    s8sym*Omega - subs(V,t,t0)]; % alternative  s8sym*subs(V,t,T) - subs(V,t,t0)

% states + constraint are SOS cones
opts.Kx = struct('sos', 0,'lin',1);
opts.Kc = struct('sos', length(prob_v_step.g));


solver_v_step = casos.sossol('S', ...
                             'mosek', ...
                             prob_v_step, ...
                             opts);

klb = casos.PS(basis(K),-inf);
kub = casos.PS(basis(K),+inf);

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

%% solve step zero
sol_s_step = solver_s_step('p',Vval, ...
                           'lbx',klb, ...
                           'ubx',kub); 

sol_v_step = solver_v_step('p',sol_s_step.x, ...
                           'lbx',Vlb, ...
                           'ubx',Vub); 


figure(1)
Omega0 = Omega;
pcontour(Omega,0,[-1 1 -1 1],'b')
hold on
pcontour(Xs,0,[-1 1 -1 1],'k--')
pcontour(subs(sol_v_step.x,t,0),0,[-1 1 -1 1],'g')
pause(0.01)


Omega_prev = subs(sol_v_step.x,t,T);
X1_tilde   = subs(sol_v_step.x,t,0);

%% setup new solver for subsequent steps
s1 = casos.PS.sym('s1',monomials([x;t],0:2));
s2 = casos.PS.sym('s2',monomials([x;t],0:2));

s3 = casos.PS.sym('s3',monomials([x;t],0:2));
s4 = casos.PS.sym('s4',monomials([x;t],0:2));

s51 = casos.PS.sym('s51',monomials([x;t],0:2));
s52 = casos.PS.sym('s52',monomials([x;t],0:2));

s61 = casos.PS.sym('s61',monomials([x;t],0:2));
s62 = casos.PS.sym('s62',monomials([x;t],0:2));

s7 = casos.PS.sym('s7',monomials(x,0:2));
s8 = casos.PS.sym('s8',monomials(x,0:2));
s9 = casos.PS.sym('s9',monomials(x,0:2));


X1_tilde_sym = casos.PS.sym('xsym',basis(V));
Omegasym     = casos.PS.sym('omegasym',basis(subs(V,t,T)));


prob_v_step = struct('x',[V;K;s1;s2;s3;s4;s51;s52;s61;s62;s7;s8;s9], ...                % decision variables
                     'f',[],...
                     'p',[X1_tilde_sym; Omegasym]); % parameter

prob_v_step.('g') = [s1;
                     s2;
                     s3;
                     s4;
                     s51;
                     s52;
                     s61;
                     s62;
                     s7;
                     s8;
                     s9
                     s1*V           - s2*hT - nabla(V,t) - nabla(V,x)*(f + gx*K); 
                     s3*V           - s4*hT - Xs;
                     K              - umin  + s51*V   - s61*hT  ; 
                     umax           - K     + s52*V   - s62*hT  ; 
                     s7*Omegasym    - subs(V,t,T);
                     s8*subs(V,t,T) - X1_tilde_sym;
                     s9*subs(V,t,T) - subs(V,t,t0)
                     ]; 

% states + constraint are SOS cones
opts = struct('sossol','mosek');
opts.Kx = struct('sos',0,'lin',length(prob_v_step.x));
opts.Kc = struct('sos', length(prob_v_step.g));
opts.verbose = 1;

tic
solver_v_step = casos.nlsossol('S1','sequential',prob_v_step,opts);
toc

Klb = casos.PS(basis(K),-inf);
Kub = casos.PS(basis(K),+inf);

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

s1lb  = casos.PS(basis(s1),-inf);
s1ub  = casos.PS(basis(s1),+inf);

s2lb  = casos.PS(basis(s2),-inf);
s2ub  = casos.PS(basis(s2),+inf);

s3lb  = casos.PS(basis(s3),-inf);
s3ub  = casos.PS(basis(s3),+inf);

s4lb  = casos.PS(basis(s4),-inf);
s4ub  = casos.PS(basis(s4),+inf);

s51lb = casos.PS(basis(s51),-inf);
s51ub = casos.PS(basis(s51),+inf);

s52lb = casos.PS(basis(s52),-inf);
s52ub = casos.PS(basis(s52),+inf);

s61lb = casos.PS(basis(s61),-inf);
s61ub = casos.PS(basis(s61),+inf);

s62lb = casos.PS(basis(s62),-inf);
s62ub = casos.PS(basis(s62),+inf);

s7lb  = casos.PS(basis(s7),-inf);
s7ub  = casos.PS(basis(s7),+inf);

s8lb  = casos.PS(basis(s8),-inf);
s8ub  = casos.PS(basis(s8),+inf);

s9lb  = casos.PS(basis(s9),-inf);
s9ub  = casos.PS(basis(s9),+inf);

lbx = [Vlb; Klb; s1lb; s2lb; s3lb; s4lb; s51lb; s52lb; s61lb; s62lb; s7lb; s8lb; s9lb];
ubx = [Vub; Kub; s1ub; s2ub; s3ub; s4ub; s51ub; s52ub; s61ub; s62ub; s7ub; s8ub; s9ub];

[~,s9_mon] = poly2basis(s9);
s9_0       = ones(length(s9_mon),1)'*s9_mon;

kstep        = 0;
NstepStorage = [cleanpoly(sol_v_step.x,1e-5)];
NstepConSet  = [cleanpoly(subs(sol_v_step.x,t,0),1e-5)];
NstepTermSet = [cleanpoly(subs(sol_v_step.x,t,T),1e-5)];

for k = 1:N

        % generate initial guess for decision variables using solution from
        % previous iteration
        if k == 1
            x0 = [sol_v_step.x;sol_s_step.x;s9_0];
        else
            x0 = sol_v_step.x;
        end

        
        sol_v_step = solver_v_step('x0',x0,'p',[X1_tilde; Omega_prev], ...
                                   'lbx',lbx, ...
                                   'ubx',ubx); 
        
        % set solutions input for next set
        Omega_prev = cleanpoly(subs(sol_v_step.x(1),t,T),1e-4);
        X1_tilde   = cleanpoly(subs(sol_v_step.x(1),t,0),1e-4);

    % store computations in array
    kstep        = [kstep; k];
    NstepStorage = [NstepStorage; cleanpoly(sol_v_step.x(1),1e-5)];
    NstepConSet  = [NstepConSet;  cleanpoly(subs(sol_v_step.x(1),t,0),1e-5)];
    NstepTermSet = [NstepTermSet; cleanpoly(subs(sol_v_step.x(1),t,T),1e-5)];
    
    % plotting of controllable sets
    figure(1)
    clf
    pcontour(Omega0,0,[-1 1 -1 1],'b')
    hold on
    pcontour(Xs,0,[-1 1 -1 1],'k--')
    pcontour(subs(sol_v_step.x(1),t,0),0,[-1 1 -1 1],'g')
    pause(0.01)
    
    disp(['Reachable Set for T: ' num2str(k*T) ' successfully computed'])
end
toc






