clc
clear
close all
profile off

%% define variables
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);
t = casos.PS('t');

x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4);
x5 = x(5);
x6 = x(6);

J = diag([8970;9230;3830]);% Casini Parameter

% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];


B = @(sigma) (1-sigma'*sigma)*eye(3)+skew(sigma)+ 2*sigma*sigma';

% dynamics
x_dot =  [-inv(J)*skew(x(1:3))*J*x(1:3) + inv(J)*u; % omega_dot
           1/4*B(x(4:6))*x(1:3)];                     % MRP kinematics

% trim point
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

[A,B] = plinearize(x_dot, x, u,x0,u0); % linearize unscaled dynamics 

% LQR controller
Q = diag([10 10 10 1 1 1]);
R = diag([0.01 0.01 0.01]);

[K,P] = lqr(A,B,Q,R);

Vval = x'*P*x;

l = Vval-0.3;
% rate constraints
omega_max =  1;
x_up  = [ omega_max omega_max omega_max]';
x_low = [-omega_max -omega_max -omega_max]';

gc   = -(x(1:3)-x_low).*(x_up-x(1:3));

% control constraints
ulow = [-1 -1 -1]';
uup  = -ulow;

umin = u-ulow;
umax = uup - u;

% time polynomial
T = 1;
hT = t*(T-t);

domain = [-1 1 -1 1]*0.05;
figure()
subplot(311)
pcontour(subs(Vval,[x1;x4;x5;x6],zeros(4,1)),1,domain)

subplot(312)
pcontour(subs(Vval,[x2;x4;x5;x6],zeros(4,1)),1,domain)

subplot(313)
pcontour(subs(Vval,[x3;x4;x5;x6],zeros(4,1)),1,domain)

domain = [-1 1 -1 1];
figure()
subplot(311)
pcontour(subs(Vval,[x4;x1;x2;x3],zeros(4,1)),1,domain)

subplot(312)
pcontour(subs(Vval,[x5;x1;x2;x3],zeros(4,1)),1,domain)

subplot(313)
pcontour(subs(Vval,[x6;x1;x2;x3],zeros(4,1)),1,domain)


% Lyapunov function candidate
V = casos.PS.sym('v',monomials([x;t],2));
k = casos.PS.sym('k',monomials([x;t],0:2),[3 1]);

% SOS multiplier
s1 = casos.PS.sym('s1',monomials([x;t],0:2),'gram');
s2 = casos.PS.sym('s2',monomials([x;t],0:2),'gram');

s3 = casos.PS.sym('s3',monomials([x;t],0:2),[3 1],'gram');
s4 = casos.PS.sym('s4',monomials([x;t],2),[3 1],'gram');

s5min = casos.PS.sym('s51',monomials([x;t],0:2),[3 1],'gram');
s6min = casos.PS.sym('s61',monomials([x;t],2),[3 1],'gram');

s5max = casos.PS.sym('s52',monomials([x;t],0:2),[3 1],'gram');
s6max = casos.PS.sym('s62',monomials([x;t],2),[3 1],'gram');

s7 = casos.PS.sym('s7',monomials(x,0:2),'gram');
s8 = casos.PS.sym('s8',monomials(x,0:2),'gram');


% level of stability
g = casos.PS.sym('g');


% options
opts               = struct('sossol','mosek');
% opts.error_on_fail = 0; 
opts.max_iter      = 20;
opts.conf_interval = [-10 0]';

opts.sossol_options.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER  = 'MSK_OPTIMIZER_INTPNT';
opts.sossol_options.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS        = 'MSK_BI_NEVER';

%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic
% profile on  -historysize 50000000000000000

V_sym = casos.PS.sym('V',basis(V));

klb = casos.PS(basis(k),-inf);
kub = casos.PS(basis(k),+inf);

% solver 1: gamma-step
sos1        = struct('x',[k;s1;s2;s3;s4;s5min;s6min;s5max;s6max;s7],...
                     'f',-g,...
                     'p',V_sym);

sos1.('g')  = [
               s1*(V_sym -g) - nabla(V_sym,x)*subs(x_dot,u,k) - s2*hT; 
               s3.*ones(3,1)*(V_sym - g) - gc - s4.*ones(3,1)*hT;
               s5min.*ones(3,1)*(V_sym - g) + subs(umin,u,k) - s6min.*ones(3,1)*hT;
               s5max.*ones(3,1)*(V_sym - g) + subs(umax,u,k) - s6max.*ones(3,1)*hT;
               s7*subs(V_sym,t,T)-l;
               ];

% states + constraint are SOS cones
opts.Kx = struct('s', 21,'l',3);
opts.Kc = struct('s', 11);
opts.verbose = 1;
S1 = casos.qcsossol('S1','bisection',sos1,opts);

tend1 = toc;
disp(['setting up solver 1 took ' num2str(tend1) ' s'])


% parameterized polynomials
s1_sym      = casos.PS.sym('s1',basis(s1));
s2_sym      = casos.PS.sym('s2',basis(s2));
s3_sym      = casos.PS.sym('s3',basis(s3(1)),[3,1]);
s4_sym      = casos.PS.sym('s41',basis(s4(1)),[3,1]);
s5min_sym   = casos.PS.sym('s51',basis(s5min(1)),[3,1]);
s5max_sym   = casos.PS.sym('s52',basis(s5max(1)),[3,1]);
s6min_sym   = casos.PS.sym('s61',basis(s6min(1)),[3,1]);
s6max_sym   = casos.PS.sym('s62',basis(s6max(1)),[3,1]);
s7_sym      = casos.PS.sym('s7',basis(s7));


Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

sos2       = struct('x',[V;s8],...
                    'p',[k;g; ...       
                         s1_sym;s2_sym;s3_sym;s4_sym; ...
                         s5min_sym;s5max_sym;s6min_sym;s6max_sym;...
                         s7_sym]);
sos2.('g') = [
               s1_sym*(V -g) - nabla(V,x)*subs(x_dot,u,k) - s2_sym*hT; 
               s3_sym.*ones(3,1)*(V - g) - gc - s4_sym.*ones(3,1)*hT;
               s5min_sym.*ones(3,1)*(V - g) + subs(umin,u,k) - s6min_sym.*ones(3,1)*hT;
               s5max_sym.*ones(3,1)*(V - g) + subs(umax,u,k) - s6max_sym.*ones(3,1)*hT;
               s7_sym*subs(V,t,T) - l;
               s8*(subs(V_sym,t,0)-g) + g - subs(V,t,0); 
             ];

opts    = struct;
opts.Kx = struct('s', 1, 'l', 1); 
opts.Kc = struct('s', 12);

opts.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
opts.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS       = 'MSK_BI_NEVER';


s2 = casos.sossol('S','mosek',sos2,opts);

tend2 = toc;

%% gamma-beta-V-iteration
disp(['setting up solver 2 took ' num2str(tend2-tend1) ' s'])
disp('==================================================================')

tic
itermax = 100; % maximum number of iterations
minIter = 5;   % minimum number of iterations (to check convergence)
disp('Start G-B-V-iteration ...')
for iter = 1:itermax

    % gamma step
    sol1 = S1('p',Vval);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            gval   = double(-sol1.f);
            s1val  = sol1.x(1);
            s3val  = sol1.x(2:4);
            s4val  = sol1.x(5:7);


            disp(['G-step feasible in ' num2str(iter) '/' num2str(itermax)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['G-step infeasible in ' num2str(iter) '/' num2str(itermax)])
                    
        otherwise
            disp(['G-step infeasible in ' num2str(iter) '/' num2str(itermax)])

    end


    % V-step
    sol3 = S3('p',[gval; s1val; s2val; s3val; s4val],...
              'lbx',Vlb,'ubx',Vub);

   switch (S3.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            Vval = sol3.x;
            disp(['V-step feasible in ' num2str(iter) '/' num2str(itermax)])
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end

    if iter <= minIter
        gvalold = gval;
    else
        if abs(gval-gvalold) <= 1e-3
            disp(['gamma is not changing anymore. Abored in iteration ' num2str(iter) '/' num2str(itermax)])
            break
        end
    end

end

tend4  = toc;
disp(['Finished G-B-V-iteration after ' num2str(tend4) ])
Vval = Vval-gval;

%% plotting
domain = [-1 1 -1 1]*0.05;
figure()
subplot(311)
pcontour(subs(Vval,[x1;x4;x5;x6],zeros(4,1)),0,domain)
hold on
pcontour(subs(p,[x1;x4;x5;x6],zeros(4,1)),0,domain,'r')

subplot(312)
pcontour(subs(Vval,[x2;x4;x5;x6],zeros(4,1)),0,domain)
hold on
pcontour(subs(p,[x2;x4;x5;x6],zeros(4,1)),0,domain,'r')
subplot(313)
pcontour(subs(Vval,[x3;x4;x5;x6],zeros(4,1)),0,domain)
hold on
pcontour(subs(p,[x3;x4;x5;x6],zeros(4,1)),0,domain,'r')


domain = [-1 1 -1 1];
figure()
subplot(311)
pcontour(subs(Vval,[x4;x1;x2;x3],zeros(4,1)),0,domain)
hold on
pcontour(subs(p,[x4;x1;x2;x3],zeros(4,1)),0,domain,'r')

subplot(312)
pcontour(subs(Vval,[x5;x1;x2;x3],zeros(4,1)),0,domain)
hold on
pcontour(subs(p,[x5;x1;x2;x3],zeros(4,1)),0,domain,'r')

subplot(313)
pcontour(subs(Vval,[x6;x1;x2;x3],zeros(4,1)),0,domain)
hold on
pcontour(subs(p,[x6;x1;x2;x3],zeros(4,1)),0,domain,'r')





