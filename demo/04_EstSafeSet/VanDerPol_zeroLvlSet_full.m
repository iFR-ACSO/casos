% Inner-approximation reachable set Van-der-Pol Oscillator

clear
close all
clc

% system states
x = casos.PS('x',2,1);
u  = casos.PS('u');
t  = casos.PS('t');
t0 = casos.PS.sym('t0');
t1 = casos.PS.sym('t1');

% system dynamics
f = [x(2);
        (1-x(1)^2)*x(2)-x(1)];

gx = [0;1];


% terminal region
P = [6.4314    0.4580
    0.4580    5.8227];

l = x'*P*x-1;              % l(x) \leq 0


% box constraint(s)
x_up  = [ 3  3];
x_low = [-3 -3];

g = [];
for k= 1:length(x_up)
    g = [g ; -(x(k)-x_low(k))*(x_up(k)-x(k) )]; % g(x) <= 0
end

% keep-out
xc = [-2;2];
lc = (x-xc)'*eye(2)*2*(x-xc)-1-1e-6; % l(x) > 0

xc = [2;-2];
lc = [lc; (x-xc)'*eye(2)*2*(x-xc)-1-1e-6]; % l(x) > 0

g = [g;-lc];

% control constraint
umin = -1;
umax = 1;

% Time horizon
T  = 2;

% time polynomial
hT = (t)*(T-t);

% Storage function and control law
V = casos.PS.sym('v',monomials([x;t],0:6));
k = casos.PS.sym('k',monomials([t;x],0:4));

% SOS multiplier

s1 = casos.PS.sym('s1',monomials([x;t],0:4),'gram'); % dissipation ineq.
s2 = casos.PS.sym('s2',monomials([x;t],4),'gram');  % dissipation ineq.

s3 = casos.PS.sym('s3',monomials([x;t],0:2),[length(g) 1],'gram');% state constraint
s4 = casos.PS.sym('s4',monomials([x;t],4),[length(g) 1],'gram');  % state constraint

s51 = casos.PS.sym('s51',monomials([x;t],0:2),'gram'); % control constraint
s52 = casos.PS.sym('s52',monomials([x;t],0:2),'gram'); % control constraint

s61 = casos.PS.sym('s61',monomials([x;t],4),'gram');   % control constraint
s62 = casos.PS.sym('s62',monomials([x;t],4),'gram');   % control constraint

s6 = casos.PS.sym('s6',monomials([x;t],4),'gram');   % control constraint 

s7 = casos.PS.sym('s7',monomials(x,0:2),'gram');   % grow constraint

%% setup solver

% solver 1: sos-multiplier and control law
tic
disp('-----------------------------------------')
disp('Building solver ... ')

% profile on
sos1 = struct('x',[k;s1;s2;s3;s4;s51;s52;s61;s62],'p',V);     % parameter

sos1.('g') = [s1*V - s2*hT - nabla(V,t) - nabla(V,x)*(f + gx*k);
              k - umin  + s51*V - s61*hT;
              umax - k   + s52*V - s62*hT;
              s3.*V.*ones(length(g),1) - s4.*hT.*ones(length(g),1) - g];

% states + constraint are SOS cones
opts    = struct;
opts.Kx = struct('s', 2+length(s3)+length(s4)+4*length(s51),'l',1);
opts.Kc = struct('s', 1+length(u)*2+length(g));

klb = casos.PS(basis(k),-inf);
kub = casos.PS(basis(k),+inf);

opts.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
opts.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';

% solver for S-step
S1 = casos.sossol('S1','mosek',sos1,opts);

% solver 2: V-step
s1sym  = casos.PS.sym('s1',basis(s1));
s2sym  = casos.PS.sym('s2',basis(s2));
s3sym  = casos.PS.sym('s3',basis(s3(1)),[length(g) 1]);
s4sym  = casos.PS.sym('s4',basis(s4(1)),[length(g) 1]);
s51sym = casos.PS.sym('s51',basis(s51));
s52sym = casos.PS.sym('s52',basis(s52));
s61sym = casos.PS.sym('s61',basis(s61));
s62sym = casos.PS.sym('s62',basis(s62));

Vsym = casos.PS.sym('V',basis(V));
ksym = casos.PS.sym('k',basis(k));


sos2 = struct('x',[V,s7]',...                                                                % dec. variable
              'p',[ksym;s1sym;s2sym;s3sym;s4sym;s51sym;s52sym;s61sym;s62sym;Vsym;t0;t1]);   % parameter

sos2.('g') = [s1sym*V   - s2sym*hT - nabla(V,t) - nabla(V,x)*(f + gx*ksym);
              ksym - umin  + s51sym*V - s61sym*hT;
              umax - ksym  + s52sym*V - s62sym*hT;
              s3sym.*V.*ones(length(g),1) - s4sym.*hT.*ones(length(g),1) - g
              subs(V,t,t1)    -  l;
              s7*subs(Vsym,t,t0)  - subs(V,t,t0)];

% states + constraint are SOS cones
opts = struct;
opts.Kx = struct('s', 1,'l',1);
opts.Kc = struct('s', 3+length(u)*2+length(g));

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

opts.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
opts.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';


% solver for V-step
S2 = casos.sossol('S2','mosek',sos2,opts);


tend1 = toc;
disp(['took ' num2str(tend1) ' s'])
% profile viewer
% profile off

%% gamma-beta-V-iteration

% initial guess storage function
Vval = x'*P*x-1;
tic

disp('-----------------------------------------')
disp('Start iteration ... ')
% profile on
for iter = 1:100

    % multiplier-step
    sol1 = S1('p',Vval,'lbx',klb,'ubx',kub);

     switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 



            disp(['S-step feasible in ' num2str(iter)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter)])
                    
        otherwise, error('Failed.')
    end


    sol2 = S2('p',[sol1.x;Vval;0;T],'lbx',Vlb,'ubx',Vub);


    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

             Vval = sol2.x(1);

             disp(['V-step feasible in ' num2str(iter)])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter)])

            tendIter = toc;
            break
                    
        otherwise, error('Failed.')
    end
  



end

tendIter = toc;


disp(['Build time is: ' num2str(tend1) ' s'])
disp(['Opt. time is: ' num2str(tendIter) ' s'])
disp('_______________________________________')
disp(['Total time is: ' num2str(tendIter+tend1) ' s'])


figure()
pcontour(subs(Vval,t,0),0,[-2 2 -2 2])
plotBoxCon([1 2],x_up,x_low)
hold on

for k= 1:length(lc)
    pcontour(lc(k),0,[-3 3 -3 3],'r--')
end
pcontour(l,0,[-2 2 -2 2],'r-')