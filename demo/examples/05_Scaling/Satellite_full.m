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

uscale =  1/(umax);

omega_max = 2*pi/180;

rateScale = 1/(omega_max-(-omega_max));

% rateScale = 180/pi;

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
Q = diag([1 1 1 1 1 1]);
R = diag([1 1 1]);

[K,P] = lqr(A,B,Q,R);

gc = [];
x_up  = [ omega_max omega_max omega_max];
x_low = [-omega_max -omega_max -omega_max];
for k= 1:length(x_up)
    gc = [gc ; -(x(k)-x_low(k))*(x_up(k)-x(k) )]; % g(x) >= 0
end

%% actual scaling
K_tilde = Du*K*Dx^(-1);  % eq.(24)

Vval = x'*(Dx^-1)'*P*(Dx^-1)*x;  % eq.(6)
Vval = cleanpoly(Vval,1e-10);

f_scaled = Dx*subs(x_dot,[x;u],[Dx^-1*x; Du^-1*(-K_tilde*x)]);  % eq.(22)


% get rid of very small terms (at least 5 magnitudes between the smallest and largest values)
f_scaled = cleanpoly(f_scaled,1e-10);

p    = Vval*2;

f = f_scaled;

u_min = u - umin; % eq. (16)
u_max = umax - u; % eq. (17)

% input is descaled
u_min = subs(u_min,u,Du^-1*u); % eq. (18)
u_max = subs(u_max,u,Du^-1*u); % eq. (19)


%% plot initial guess
lvlSet = 1;

domain = [-1 1 -1 1]*0.05;
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
s5 = casos.PS.sym('s3',monomials(x,0:2),[length(x_up) 1],'gram');


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
opts.conf_interval = [-10 0]';


%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic

V_sym = casos.PS.sym('V',basis(V));

% solver 1: gamma-step
sos1        = struct('x',[s1;s3;s4;s5], ...
                     'f',-g, ...
                     'p',V_sym);

sos1.('g')  = [s1*(V_sym - g) - nabla(V_sym , x)*f - l;
               subs(u_min,u,-K_tilde*x) + s3*(V_sym - g);
               subs(u_max,u,-K_tilde*x) + s4*(V_sym - g)
               s5*(V_sym-g) - subs(gc,x,Dx^-1*x)];   % eq.(37)

% states + constraint are SOS cones
opts.Kx = struct('s', 1+length(s4)+length(s3)+length(s5));
opts.Kc = struct('s', 1+length(s4)+length(s3)+length(s5));

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
opts.error_on_fail  = 0; 
opts.max_iter       = 20;
opts.verbose        = 1;
opts.conf_interval  = [-10 0]';

S2 = casos.qcsossol('S2','bisection',sos2,opts);

tend2 = toc;
disp(['setting up solver 2 took ' num2str(tend2-tend1) ' s'])

% solver 3: V-step
s1_sym = casos.PS.sym('s1',basis(s1));
s2_sym = casos.PS.sym('s2',basis(s2));
s3_sym = casos.PS.sym('s3',basis(s3(1)),[length(s3) 1]);
s4_sym = casos.PS.sym('s4',basis(s4(1)),[length(s4) 1]);
s5_sym = casos.PS.sym('s5',basis(s5(1)),[length(s5) 1]);

% s1 = casos.PS.sym()
Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

sos3       = struct('x',V, ...
                    'p',[s1_sym;s3_sym; s4_sym; s5_sym; s2_sym; g; b]); % order intentionally

sos3.('g') = [V - l;
              s1_sym*(V - g) - nabla(V,x)*f - l
              s2_sym*(p - b) + g - V;
              subs(u_min,u,-K_tilde*x) + s3_sym*(V - g);
              subs(u_max,u,-K_tilde*x) + s4_sym*(V - g)
              s5_sym*(V - g) - subs(gc,x,Dx^-1*x)];  % eq.(37)

opts    = struct;
opts.Kx = struct('s', 0, 'l', 1); 
opts.Kc = struct('s', 3+length(s3) +length(s4) +length(s5));

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
    sol3 = S3('p',[sol1.x;sol2.x;gval;bval],'lbx',Vlb,'ubx',Vub);

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




Vval = subs(Vval-gval,x,Dx*x);  % eq.(13)
% Vval = subs(Vval,x(1:3),x(1:3).*[180/pi;180/pi;180/pi]);
p = subs(p-bval,x,Dx*x);
% p = subs(p,x(1:3),x(1:3).*[180/pi;180/pi;180/pi]);

%% plotting
domain = [-1 1 -1 1]*1.5*omega_max;
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





