% Inner-approximation reachable set Van-der-Pol Oscillator

clear
close all
clc

%% define variables
x = casos.PS('x',4,1);
u = casos.PS('u',1,1);

t  = casos.PS('t');
t0 = casos.PS.sym('t0');
t1 = casos.PS.sym('t1');

x1 = x(1);
x2 = x(2); 
x3 = x(3); 
x4 = x(4); 

% scaled 4-state system
f2 = - 10.6560*x1^3 + 11.5309*x1^2*x3 + 7.8850*x1*x3^2 + 0.7972*x2^2*x3 ...
  + 0.8408*x2*x3*x4 + 21.0492*x3^3 + 0.4204*x3*x4^2 + 66.5225*x1 - 24.5110*x3;

f4 = 10.9955*x1^3 - 48.9151*x1^2*x3 - 6.4044*x1*x3^2 - 2.3955*x2^2*x3 ...
  - 1.5943*x2*x3*x4 - 51.9088*x3^3 - 0.7971*x3*x4^2 - 68.6419*x1 + 103.9783*x3;

f = [x2; f2; x4; f4];

g2 = -10.0959*x3^2 + 44.2521;
g4 =  37.8015*x3^2 - 83.9120;

gx = [0; g2; 0; g4];


% target function
l = x'*blkdiag(1/0.1^2, 1/0.35^2, 1/0.1^2, 1/0.35^2)*x - 1;

% l =x'*x-1;
Vval = 182.2102*x1^2 + 67.8981*x1*x2 + 314.3265*x1*x3 + 37.2705*x1*x4 + ...
    6.4123*x2^2 + 59.0528*x2*x3 + 7.0314*x2*x4 + 138.9343*x3^2 + ...
    32.5944*x3*x4 + 1.9521*x4^2;

% control constraint
umin = -1;
umax = 1;

% Time horizon
T  = 4;

figure()
pcontour((subs(Vval,[t;x2;x4],zeros(3,1))),0,[-0.5 0.5 -0.5 0.5]*2)
hold on
pcontour(subs(l,[t;x3;x4],zeros(3,1)),0,[-1 1 -1 1],'r-')

% time polynomial
hT = (t)*(T-t);

% Storage function and control law
V = casos.PS.sym('v',monomials([x;t],0:4));
k = casos.PS.sym('k',monomials([t;x],0:4));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials([x;t],0:4)); % dissipation ineq.
s2 = casos.PS.sym('s2',monomials([x;t],0:4));  % dissipation ineq.

s3 = casos.PS.sym('s3',monomials([x;t],0:4));% state constraint
s4 = casos.PS.sym('s4',monomials([x;t],0:4));  % state constraint

s51 = casos.PS.sym('s51',monomials([x;t],0:4)); % control constraint
s52 = casos.PS.sym('s52',monomials([x;t],0:4)); % control constraint

s61 = casos.PS.sym('s61',monomials([x;t],0:4));   % control constraint
s62 = casos.PS.sym('s62',monomials([x;t],0:4));   % control constraint

s7 = casos.PS.sym('s7',monomials(x,0:4));   % grow constraint
s8 = casos.PS.sym('s8',monomials(x,0:4));   % grow constraint


%% setup solver

% solver 1: sos-multiplier and control law
tic
disp('-----------------------------------------')
disp('Building solver ... ')
gamma = casos.PS.sym('g');

% profile on
sos1 = struct('x',[k,s1,s2,s51,s52,s61,s62,s8]','p',[V;gamma]);     % parameter

sos1.('g') = [s1;s2;s51;s52;s61;s62;s8-0.0001;
              s1*(V-gamma) - s2*hT - nabla(V,t) - nabla(V,x)*(f + gx*k);
              subs(V,t,T)-gamma    -  s8*l;
              k - umin  + s51*(V-gamma) - s61*hT;
              umax - k   + s52*(V-gamma) - s62*hT];

% states + constraint are SOS cones
opts    = struct;
opts.error_on_fail = 0;
opts.Kx = struct('s', 0,'l',8);
opts.Kc = struct('s', length(sos1.g));

klb = casos.PS(basis(k),-inf);
kub = casos.PS(basis(k),+inf);


% solver for S-step
S1 = casos.sossol('S1','mosek',sos1,opts);

% solver 2: V-step
s1sym  = casos.PS.sym('s1',basis(s1));
s2sym  = casos.PS.sym('s2',basis(s2));
s3sym  = casos.PS.sym('s3',basis(s3));
s4sym  = casos.PS.sym('s4',basis(s4));
s51sym = casos.PS.sym('s51',basis(s51));
s52sym = casos.PS.sym('s52',basis(s52));
s61sym = casos.PS.sym('s61',basis(s61));
s62sym = casos.PS.sym('s62',basis(s62));
s8sym = casos.PS.sym('s8',basis(s8));

Vsym = casos.PS.sym('V',basis(V));
ksym = casos.PS.sym('k',basis(k));


sos2 = struct('x',[V;s7;s8],...
              'p',[gamma,s1sym,s2sym,s51sym,s52sym,s61sym,s62sym,ksym,Vsym,t0,t1]');   % parameter

sos2.('g') = [s7;s8-0.0001;
              s1sym*(V-gamma)   - s2sym*hT - nabla(V,t) - nabla(V,x)*(f + gx*ksym);
              ksym - umin  + s51sym*(V-gamma) - s61sym*hT;
              umax - ksym  + s52sym*(V-gamma) - s62sym*hT;
              (subs(V,t,T) - gamma)    -  s8*l;
              s7*(subs(Vsym,t,0) - gamma) - subs(V,t,0) + gamma];

% states + constraint are SOS cones
opts = struct;
opts.error_on_fail = 0;
opts.Kx = struct('s', 0,'l',3);
opts.Kc = struct('s', length(sos2.g));

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

% solver for V-step
S2 = casos.sossol('S2','mosek',sos2,opts);


tend1 = toc;
disp(['took ' num2str(tend1) ' s'])


%% gamma-beta-V-iteration

% initial guess storage function
tic

disp('-----------------------------------------')
disp('Start iteration ... ')
gamma_use = [];
for iter = 1:10

    % Initialize the bisection
    gamma_ub = 1;
    gamma_lb = 0;
    num_exp = 0;
    go = 1;
    while (num_exp <= 12 || go)

        num_exp = num_exp + 1;
        if num_exp >= 20
            break
        end

    gamma_try = (gamma_ub + gamma_lb)/2 

    sol1 = S1('p',[Vval;gamma_try],'lbx',klb,'ubx',kub);

     switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            kval   = sol1.x(1);
            s1val  = sol1.x(2);
            s2val  = sol1.x(3);
            s51val = sol1.x(4);
            s52val = sol1.x(5);
            s61val = sol1.x(6);
            s62val = sol1.x(7);
            % s8val  = sol1.x(8);
            gamma_lb = gamma_try;
            go = 0;
            gamma_use = gamma_try;

            disp(['S-step feasible in ' num2str(iter)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter)])
                    gamma_ub = gamma_try; 
         otherwise
              gamma_ub = gamma_try;
             % break
     end
    end

    if isempty(gamma_use)
        break
    end

    sol2 = S2('p',[gamma_use,s1val,s2val,s51val,s52val,s61val,s62val,kval,Vval,0,T],'lbx',Vlb,'ubx',Vub);


    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

             Vval = sol2.x(1);

             disp(['V-step feasible in ' num2str(iter)])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter)])

                break
        otherwise
            disp(['V-step infeasible in ' num2str(iter)])
            break
    end
  


end

tendIter = toc;

%%
figure()
pcontour((subs(Vval,[t;x2;x4],zeros(3,1))),gamma_use,[-0.5 0.5 -0.5 0.5])
hold on
pcontour((subs(Vval,[t;x2;x4],[T;zeros(2,1)])),gamma_use,[-0.5 0.5 -0.5 0.5],'g--')
pcontour(subs(l,[t;x2;x4],zeros(3,1)),0,[-0.5 0.5 -0.5 0.5],'r-')
legend('V(0,x)','V(T,x)','l(x)')
