clc
clear
close all
profile off

%% define variables
x = casos.PS('x',6,1);
u = casos.PS('u',3,1);

x1 = x(1); 
x2 = x(2); 
x3 = x(3); 
x4 = x(4);
x5 = x(5);
x6 = x(6);

J = diag([8970;9230;3830]);% Casini Parameter

umax =  0.16; % Nm or 160 mNm
umin = -umax;

uscale = 1/(umax-umin);

omega_max = 2*pi/180;

rateScale = 1/(omega_max-(-omega_max));

Dx = diag([rateScale,rateScale,rateScale,1,1,1]);
Du = diag([uscale,uscale,uscale]);


% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];


B = @(sigma) (1-sigma'*sigma)*eye(3)+skew(sigma)+ 2*sigma*sigma';

% dynamics
x_dot =  [-inv(J)*skew(x(1:3))*J*x(1:3) + inv(J)*u; % omega_dot
           1/4*B(x(4:6))*x(1:3)];                   % MRP kinematics

% trim point
x0    = [0 0 0 0 0 0]';
u0    = [0,0,0]';

[A,B] = plinearize(x_dot, x, u,x0,u0); % linearize unscaled dynamics 

% LQR controller
Q = diag([0.1 0.1 0.1 0.1 0.1 0.1]);
R = diag([1 1 1]);

[K,P] = lqr(A,B,Q,R);


%% actual scaling
K_tilde = Du*K*Dx^(-1);  % eq.(24)

Vval = x'*(Dx^-1)'*P*(Dx^-1)*x;  % eq.(6)
Vval = cleanpoly(Vval,1e-10);

f_scaled = Dx*subs(x_dot,[x;u],[Dx^-1*x; Du^-1*(-K_tilde*x)]); % eq. (22)

u_min = u - umin; 
u_max = umax - u;


u_min = subs(u_min,u,Du^-1*u); % eq. (18)
u_max = subs(u_max,u,Du^-1*u); % eq. (19)


% get rid of very small terms (at least 5 magnitudes)
f_scaled = cleanpoly(f_scaled,1e-10);

p    = Vval*2;

f = f_scaled;

lvlSet = 1;

domain = [-1 1 -1 1];
figure()
subplot(311)
pcontour(subs(Vval,[x1;x4;x5;x6],zeros(4,1)),lvlSet,domain)
hold on
pcontour(subs(p,[x1;x4;x5;x6],zeros(4,1)),lvlSet,domain,'r')

subplot(312)
pcontour(subs(Vval,[x2;x4;x5;x6],zeros(4,1)),lvlSet,domain)
hold on
pcontour(subs(p,[x2;x4;x5;x6],zeros(4,1)),lvlSet,domain,'r')
subplot(313)
pcontour(subs(Vval,[x3;x4;x5;x6],zeros(4,1)),lvlSet,domain)
hold on
pcontour(subs(p,[x3;x4;x5;x6],zeros(4,1)),lvlSet,domain,'r')

figure()
subplot(311)
pcontour(subs(Vval,[x4;x1;x2;x3],zeros(4,1)),lvlSet,domain)
hold on
pcontour(subs(p,[x4;x1;x2;x3],zeros(4,1)),lvlSet,domain,'r')

subplot(312)
pcontour(subs(Vval,[x5;x1;x2;x3],zeros(4,1)),lvlSet,domain)
hold on
pcontour(subs(p,[x5;x1;x2;x3],zeros(4,1)),lvlSet,domain,'r')

subplot(313)
pcontour(subs(Vval,[x6;x1;x2;x3],zeros(4,1)),lvlSet,domain)
hold on
pcontour(subs(p,[x6;x1;x2;x3],zeros(4,1)),lvlSet,domain,'r')


% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,2),'gram');
s2 = casos.PS.sym('s2',monomials(x,2),'gram');

s3 = casos.PS.sym('s3',monomials(x,0:2),[length(u) 1],'gram');
s4 = casos.PS.sym('s4',monomials(x,0:2),[length(u) 1],'gram');

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
g = casos.PS.sym('g');
b = casos.PS.sym('b');

% options
opts               = struct('sossol','mosek');
opts.error_on_fail = 0; 
opts.max_iter      = 20;
opts.verbose       = 1;
opts.conf_interval = [-1000 0]';


%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic

V_sym = casos.PS.sym('V',basis(V));

% solver 1: gamma-step
sos1        = struct('x',[s1;s3;s4], ...
                     'f',-g, ...
                     'p',V_sym);

sos1.('g')  = [s1*(V_sym - g) - nabla(V_sym , x)*f - l;
               u_min + s3.ones(length(u),1).*(V_sym - g);
               u_max + s4.ones(length(u),1).*(V_sym - g)];

% states + constraint are SOS cones
opts.Kx = struct('s', 1+length(s3)+length(s4));
opts.Kc = struct('s', 1+length(s3)+length(s4));

S1 = casos.qcsossol('S1','bisection',sos1,opts);

tend1 = toc;
disp(['setting up solver 1 took ' num2str(tend1) ' s'])


% solver 2: beta-step
sos2        = struct('x',s2, ...
                     'f',-b, ...
                     'p',[V_sym;g]);

sos2.('g')  = s2*(p - b) + g - V_sym;

% states + constraint are SOS cones
opts.Kx             = struct('s', 1);
opts.Kc             = struct('s', 1);
opts.error_on_fail  = 1; 
opts.max_iter       = 20;
opts.verbose        = 1;
opts.conf_interval  = [-1000 0]';

S2 = casos.qcsossol('S2','bisection',sos2,opts);

tend2 = toc;
disp(['setting up solver 2 took ' num2str(tend2-tend1) ' s'])

% solver 3: V-step
s1_sym = casos.PS.sym('s1',basis(s1));
s2_sym = casos.PS.sym('s2',basis(s2));
s3_sym = casos.PS.sym('s3',basis(s3(1)),[length(s3) 1]);
s4_sym = casos.PS.sym('s4',basis(s2(1)),[length(s4) 1]);

% s1 = casos.PS.sym()
Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

sos3       = struct('x',V, ...
                    'p',[s1_sym; s2_sym; s3_sym; s4_sym; g; b]);

sos3.('g') = [V - l;
              s1_sym*(V - g) - nabla(V,x)*f - l
              s2_sym*(p - b) + g - V;
              u_min + s3_sym.ones(length(u),1).*(V - g);
              u_max + s4_sym.ones(length(u),1).*(V - g)];

opts    = struct;
opts.Kx = struct('s', 0, 'l', 1); 
opts.Kc = struct('s', 3+length(s3)+length(s4));

S3 = casos.sossol('S','mosek',sos3,opts);

tend3 = toc;

%% gamma-beta-V-iteration
disp(['setting up solver 3 took ' num2str(tend3-tend2) ' s'])
disp('==================================================================')

tic
itermax = 10;
disp('Start G-B-V-iteration ...')
for iter = 1:itermax

    % gamma step
    sol1 = S1('p',Vval);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            gval = double(-sol1.f);
            % s1val = sol1.x;
            disp(['G-step feasible in ' num2str(iter) '/' num2str(itermax)])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['G-step infeasible in ' num2str(iter) '/' num2str(itermax)])
                    return
        otherwise, error('Failed.')
    end

    % beta step
    sol2 = S2('p',[Vval;gval]);
   switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            bval = -sol2.f;
            s2val = sol2.x;

            disp(['B-step feasible in ' num2str(iter) '/' num2str(itermax) ])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['B-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end

    % V-step
    sol3 = S3('p',[sol1.x;gval;bval],'lbx',Vlb,'ubx',Vub);

   switch (S3.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            Vval = sol3.x;

            disp(['V-step feasible in ' num2str(iter) '/' num2str(itermax)])
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end
end

tend4  = toc;
disp(['Finished G-B-V-iteration after ' num2str(tend4) ])


Vval = subs(Vval-gval,x,Dx*x); % eq.(13)
p = subs(p-bval,x,Dx*x);

%% plotting
domain = [-1 1 -1 1]*0.5;
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





