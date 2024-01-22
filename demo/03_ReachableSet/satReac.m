% Inner-approximation reachable set Van-der-Pol Oscillator

clear
clc

% indeterminates
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);
t = casos.PS('t');

x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4);
x5 = x(5);
x6 = x(6);

nx  =  length(x);
nu  = length(u);


J = diag([8970;9230;3830]);% Casini Parameter

% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];


B = @(sigma) (1-sigma'*sigma)*eye(3)+skew(sigma)+ 2*sigma*sigma';

% dynamics
x_dot =  [-inv(J)*skew(x(1:3))*J*x(1:3) + inv(J)*u; % omega_dot
           1/4*B(x(4:6))*x(1:3)];  

% rate constraints
omega_max =  1;
x_up  = [ omega_max omega_max omega_max]';
x_low = [-omega_max -omega_max -omega_max]';

gc   = -(x(1:3)-x_low).*(x_up-x(1:3));

% control constraints
umin = [-1 -1 -1]';
umax  = -umin;


ng = length(gc);

l = -0.269318 + 15838.7*x(1)^2 + 16403.9*x(2)^2 + 5433.14*x(3)^2 + 1480.93*x(1)...
  *x(4) + 73.3*x(4)^2 + 1512.43*x(2)*x(5) + 74.0522*x(5)^2 + 765.567*x(3)*x(6) ...
  + 53.4006*x(6)^2;

% Time horizon
T     = 1;

% time polynomial
hT = (t)*(T-t);


% Storage function and control law
V = casos.PS.sym('v',monomials([x;t],0:2));
k = casos.PS.sym('k',monomials([x;t],0:2),[nu 1]);

% SOS multiplier

s1 = casos.PS.sym('s1',monomials([x;t],0:4),'gram'); % dissipation ineq.
s2 = casos.PS.sym('s2',monomials([x;t],4),'gram');  % dissipation ineq.

s3 = casos.PS.sym('s3',monomials([x;t],0:2),[ng 1],'gram');  % state constraint
s4 = casos.PS.sym('s4',monomials([x;t],0:2),[ng 1],'gram');  % state constraint

s51 = casos.PS.sym('s51',monomials([x;t],0:2),[nu 1],'gram'); % control constraint
s52 = casos.PS.sym('s52',monomials([x;t],0:2),[nu 1],'gram'); % control constraint

s61 = casos.PS.sym('s61',monomials([x;t],2),[nu 1],'gram');   % control constraint
s62 = casos.PS.sym('s62',monomials([x;t],2),[nu 1],'gram');   % control constraint

s7 = casos.PS.sym('s7',monomials(x,0:2),'gram');   % terminal set containment
s8 = casos.PS.sym('s8',monomials(x,0:2),'gram');   % grow constraint

%% setup solver

% solver 1: sos-multiplier and control law
tic
disp('-----------------------------------------')
disp('Building solver ... ')

% profile on
g = casos.PS.sym('g');

sos1 = struct('x',[k;s1;s2;s3;s4;s51;s52;s61;s62;s7]... % decsision variables
             ,'f',-g,...   % cost
              'p',V');     % parameter

sos1.('g') = [s1*(V-g) - s2*hT - nabla(V,t) - nabla(V,x)*subs(x_dot,u,k);
              k - umin  + s51.*ones(nu,1)*(V-g) - s61.*ones(nu,1)*hT;
              umax - k   + s52.*ones(nu,1)*(V-g) - s62.*ones(nu,1)*hT;
              s7*(subs(V,t,T)-g)    -  l;
              s3.*ones(ng,1)*(V-g) - s4.*ones(ng,1)*hT - gc];

% states + constraint are SOS cones
opts = struct('sossol','mosek');
opts.Kx = struct('s', 2+2*ng+4*nu+1,'l',nu);
opts.Kc = struct('s', 2+ng+2*nu);
opts.error_on_fail = 0; 
opts.verbose = 1;
opts.conf_interval = [-10 0]';

klb = casos.PS(basis(k(1)),-inf)*ones(nu,1);
kub = casos.PS(basis(k(1)),+inf)*ones(nu,1);

opts.sossol_options.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
opts.sossol_options.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';

% solver for S-step
S1 = casos.qcsossol('S1','bisection',sos1,opts);

tend1 = toc;
disp(['took ' num2str(tend1) ' s'])


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
ksym = casos.PS.sym('k',basis(k(1)),[nu 1]);


sos2 = struct('x',[V,s8]',...                                                                        % dec. variable
              'p',[ksym;s1sym;s2sym;s3sym;s4sym;s51sym;s52sym;s61sym;s62sym;s7sym;Vsym;g]);   % parameter

sos2.('g') = [s1sym*(V-g)  - s2sym*hT     - nabla(V,t) - nabla(V,x)*subs(x_dot,u,ksym);
              ksym - umin  + s51sym.*ones(nu,1)*(V-g) - s61sym.*ones(nu,1)*hT;
              umax - ksym  + s52sym.*ones(nu,1)*(V-g) - s62sym.*ones(nu,1)*hT;
              s3sym.*ones(ng,1)*(V-g)  - s4sym.*ones(ng,1)*hT - gc;
              s7sym*(subs(V,t,T)-g)    -  l;
              s8*(subs(Vsym,t,0)-g)    + g - subs(V,t,0)];

% states + constraint are SOS cones
opts = struct;
opts.Kx = struct('s', 1,'l',1);
opts.Kc = struct('s', 3+ng+2*nu);

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

opts.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
opts.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';


% solver for V-step
S2 = casos.sossol('S2','mosek',sos2,opts);


tend2 = toc;
disp(['took ' num2str(tend2-tend1) ' s'])


%% gamma-beta-V-iteration

% initial guess storage function
Vval = l;

tic

disp('-----------------------------------------')
disp('Start iteration ... ')

for iter = 1:100

    % multiplier-step
    sol1 = S1('p',Vval','lbx',klb,'ubx',kub);

     switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            % kval   = sol1.x(1);
            % s1val  = sol1.x(2);
            % s2val  = sol1.x(3);
            % s3val  = sol1.x(4);
            % s4val  = sol1.x(5);
            % s51val = sol1.x(6);
            % s52val = sol1.x(7);
            % s61val = sol1.x(8);
            % s62val = sol1.x(9);
            % s7val = sol1.x(10);
            gval  = double(-sol1.f);

            disp(['S-step feasible in ' num2str(iter)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter)])
                    
        otherwise, error('Failed.')
    end

    % V-step
    sol2 = S2('p',[sol1.x ;Vval;gval],'lbx',Vlb,'ubx',Vub);


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

tendIter = toc;


%% plotting
figure()
pcontour(subs(Vval,t,0)-gval,0,[-2 2 -2 2])
hold on
pcontour((subs(Vval,t,1)),0,[-2 2 -2 2],'k--')
pcontour(gc(1),0,[-2 2 -2 2],'k-')

pcontour(l,0,[-2 2 -2 2],'r-')