% Inner-approximation reachable set Van-der-Pol Oscillator

clear
close all
clc
profile off


% indeterminates
x = casos.PS('x',2,1);
w = casos.PS('w',1,1);
t = casos.PS('t');


nx  =  length(x);
nu  = 1;

% system dynamics
f = [x(2);
        (1-x(1)^2)*x(2)-x(1)];

gx = [0;1];

gw = [0;1];


% terminal region
P = [6.4314    0.4580
     0.4580    5.8227];

l = x'*P*x-1;              % l(x) \leq 0

% state constraint
gc = 3*x(2)^2 + x(1)^2-1;  % g(x) \leq 0


ng = length(gc);


% control constraint
umin = -1;
umax =  1;

% Time horizon
T     = 1;

% time polynomial
hT = (t)*(T-t);


% Storage function and control law
V = casos.PS.sym('v',monomials([x;t],0:4));
k = casos.PS.sym('k',monomials([x;t],0:4));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials([x;t;w],0:2),'gram'); % dissipation ineq.
s2 = casos.PS.sym('s2',monomials([x;t;w],2),'gram');  % dissipation ineq.

s3 = casos.PS.sym('s3',monomials([x;t;w],0:2),[ng 1],'gram');  % state constraint
s4 = casos.PS.sym('s4',monomials([x;t;w],2),[ng 1],'gram');  % state constraint

s51 = casos.PS.sym('s51',monomials([x;t;w],0:2),[nu 1],'gram'); % control constraint
s52 = casos.PS.sym('s52',monomials([x;t;w],0:2),[nu 1],'gram'); % control constraint

s61 = casos.PS.sym('s61',monomials([x;t;w],2),[nu 1],'gram');   % control constraint
s62 = casos.PS.sym('s62',monomials([x;t;w],2),[nu 1],'gram');   % control constraint

s7 = casos.PS.sym('s7',monomials(x,0:2),'gram');   % terminal set containment
s8 = casos.PS.sym('s8',monomials(x,0:2),'gram');   % grow constraint


s9  = casos.PS.sym('s9',monomials([x;t;w],2),'gram');   % terminal set containment
s10 = casos.PS.sym('s10',monomials([x;t;w],2),'gram');   % grow constraint
s11 = casos.PS.sym('s11',monomials([x;t;w],2),'gram');   % grow constraint


%% setup solver
profile on -historysize 500000000000 
% solver 1: sos-multiplier and control law
tic
disp('-----------------------------------------')
disp('Building solver ... ')

% profile on
g = casos.PS.sym('g');

R = 0.1;
q = 1;

sos1 = struct('x',[k;s1;s2;s3;s4;s51;s52;s61;s62;s7;s9;s10;s11]... % decsision variables
             ,'f',-g,...   % cost
              'p',V);     % parameter

paraLght = length(sos1.x)-size(s9,1)-size(s10,1)-size(s11,1);

sos1.('g') = [s1*(V - g - R^2*q) - s2*hT - nabla(V,t) - nabla(V,x)*(f + gx*k + gw*w)  + w'*w  +  s9*(w'*w - 2*R^2);
              k - umin   + s51.*ones(nu,1)*(V - g - R^2*q) - s61.*ones(nu,1)*hT        + s10*(w'*w - 2*R^2);
              umax - k   + s52.*ones(nu,1)*(V - g - R^2*q) - s62.*ones(nu,1)*hT        + s11*(w'*w - 2*R^2);
              s7*(subs(V,t,T) - g - R^2*q)    -  l;
              s3.*ones(ng,1)*(V - g - R^2*q) - s4.*ones(ng,1)*hT - gc];

% states + constraint are SOS cones
opts = struct('sossol','mosek');
opts.Kx = struct('s', 2+2*ng+4*nu+1+3,'l',1);
opts.Kc = struct('s', 2+ng+2*nu);
opts.error_on_fail = 0; 
opts.verbose = 1;
opts.conf_interval = [-10 0]';

klb = casos.PS(basis(k),-inf);
kub = casos.PS(basis(k),+inf);

opts.sossol_options.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
opts.sossol_options.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';

% solver for S-step
S1 = casos.qcsossol('S1','bisection',sos1,opts);

% solver 2: V-step
s1sym  = casos.PS.sym('s1',basis(s1));
s2sym  = casos.PS.sym('s2',basis(s2));
s3sym  = casos.PS.sym('s3',basis(s3(1)),[ng 1]);
s4sym  = casos.PS.sym('s4',basis(s4(1)),[ng 1]);
s51sym = casos.PS.sym('s51',basis(s51(1)),[nu 1]);
s52sym = casos.PS.sym('s52',basis(s52(1)),[nu 1]);
s61sym = casos.PS.sym('s61',basis(s61(1)),[nu 1]);
s62sym = casos.PS.sym('s62',basis(s62(1)),[nu 1]);
s7sym  = casos.PS.sym('s7',basis(s7));


Vsym = casos.PS.sym('V',basis(V));
ksym = casos.PS.sym('k',basis(k));


sos2 = struct('x',[V;s8;s9;s10;s11],...                                                                        % dec. variable
              'p',[ksym;s1sym;s2sym;s3sym;s4sym;s51sym;s52sym;s61sym;s62sym;s7sym;Vsym;g]);   % parameter

sos2.('g') = [s1sym*(V - g - R^2*q)  - s2sym*hT     - nabla(V-g,t) - nabla(V-g,x)*(f + gx*ksym+ gw*w) + w'*w   +  s9*(w'*w - 2*R^2) ;
              ksym - umin  + s51sym.*ones(nu,1)*(V - g - R^2*q) - s61sym.*ones(nu,1)*hT              + s10*(w'*w - 2*R^2);
              umax - ksym  + s52sym.*ones(nu,1)*(V - g - R^2*q) - s62sym.*ones(nu,1)*hT              + s11*(w'*w - 2*R^2);
              s3sym.*ones(ng,1)*(V - g - R^2*q)  - s4sym.*ones(ng,1)*hT - gc;
              s7sym*(subs(V,t,T) - g - R^2*q)    -  l;
              s8*(subs(Vsym,t,0) - g - R^2*q)    - (subs(V,t,0) - g - R^2*q)];

% states + constraint are SOS cones
opts = struct;
opts.Kx = struct('s', 4,'l',1);
opts.Kc = struct('s', 3+ng+2*nu);

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

opts.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
opts.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS       = 'MSK_BI_NEVER';


% solver for V-step
S2 = casos.sossol('S2','mosek',sos2,opts);


tend1 = toc;
disp(['took ' num2str(tend1) ' s'])




%% gamma-beta-V-iteration

% initial guess storage function
Vval = x'*P*x;

tic

disp('-----------------------------------------')
disp('Start iteration ... ')

for iter = 1:10

    % multiplier-step
    sol1 = S1('p',Vval','lbx',klb,'ubx',kub);

     switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            gval  = double(-sol1.f);

            disp(['S-step feasible in ' num2str(iter)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter)])
                    
        otherwise, error('Failed.')
    end

    % V-step
    sol2 = S2('p',[sol1.x(1:paraLght); Vval;gval],'lbx',Vlb,'ubx',Vub);


    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

             Vval = sol2.x(1);

             disp(['V-step feasible in ' num2str(iter)])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter)])


            break
                    
        otherwise, error('Failed.')
    end
  



end
tBuild =tend1;
tendIter = toc;

tTotal = tBuild+tendIter;

disp(['Build time is ' num2str(tBuild) ' s'])
disp(['Itertaion time is ' num2str(tendIter) ' s'])
disp('______________________________________________')
disp(['Total time is ' num2str(tTotal) ' s'])


profile viewer

%% plotting
figure()
pcontour(subs(Vval,t,0)-gval,0,[-2 2 -2 2],'b')
hold on
pcontour(subs(Vval,t,1)-gval,0,[-2 2 -2 2],'b--')
pcontour(gc(1),0,[-2 2 -2 2],'k-')

pcontour(l,0,[-2 2 -2 2],'r-')