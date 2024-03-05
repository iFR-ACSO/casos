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

load TermSet.mat

% time
t  = casos.PS('t');
t0 = casos.PS.sym('t0');
t1 = casos.PS.sym('t1');


g  = casos.PS.sym('g');

% reachability storage function
V = casos.PS.sym('v',monomials([x;t],0:4));
k = casos.PS.sym('k',monomials([x;t],0:3),[length(u) 1]);

% SOS multiplier
s1 = casos.PS.sym('s1',monomials([x;t],0:4));
s2 = casos.PS.sym('s2',monomials([x;t],0:4));

s3 = casos.PS.sym('s3',monomials([x;t],0:4),[3 1]);
s4 = casos.PS.sym('s4',monomials([x;t],0:3),[3 1]);

s51 = casos.PS.sym('s51',monomials([x;t],0:3),[length(u) 1]);
s52 = casos.PS.sym('s52',monomials([x;t],0:3),[length(u) 1]);

s61 = casos.PS.sym('s61',monomials([x;t],0:3),[length(u) 1]);
s62 = casos.PS.sym('s62',monomials([x;t],0:3),[length(u) 1]);

s7 = casos.PS.sym('s7',monomials(x,0:4));   % terminal-set inclusion
s8 = casos.PS.sym('s8',monomials(x,0:4));   % grow constraint

% Time horizon
T  = 0.1;

% time polynomial
hT = (t)*(T-t);


%% define satellite
J = diag([8970;9230;3830]);% Casini Parameter

umax =  0.16; % Nm or 160 mNm
umin = -umax;

uscale =  1/(umax);

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

x_up  = [ omega_max  omega_max  omega_max]';
x_low = [-omega_max -omega_max -omega_max]';

gc = -(x(1:3)-x_low).*(x_up-x(1:3)) ;


u_min = u - umin; % eq. (16)
u_max = umax - u; % eq. (17)

load TermSet.mat
Vval = lscaled;
l    =  Vval;

% x_dot = Dx*subs(x_dot,x,Dx^-1*x);

%% setup solver

profile on -historysize 50000000000000
disp('==================================================================')
disp('Setting up solver ... ')
tic

V_sym = casos.PS.sym('V',basis(V));

sos1 = struct('x',[k;s1;s3;s51;s52;s2;s4;s61;s62],...
              'p',V_sym);     

sos1.('g') = [s1;s3;s51;s52;s2;s4;s61;s62;
              s1*(V_sym) - s2*hT - nabla((V_sym),t) - nabla((V_sym),x)*subs(x_dot,u,k);
              s3.*(V_sym) - s4.*hT  - subs(gc,x,Dx^-1*x);
              subs(u_min,u,k)     + s51.*(V_sym) - s61.*hT;
              subs(u_max,u,k)     + s52.*(V_sym) - s62.*hT];

% states + constraint are SOS cones
opts    = struct;
opts.Kx = struct('s', 0, ... % size(sos1.x((length(u)+1):end),1)
                 'l',size(sos1.x,1));
opts.Kc = struct('s', length(sos1.g));

S1 = casos.sossol('S1','mosek',sos1,opts);

tend1 = toc;
disp(['setting up solver 1 took ' num2str(tend1) ' s'])

profile viewer

s1_sym = casos.PS.sym('s1',basis(s1));
s2_sym = casos.PS.sym('s2',basis(s1));
s3_sym = casos.PS.sym('s3',basis(s3(1)),[3 1]);
s4_sym = casos.PS.sym('s4',basis(s4(1)),[3 1]);

s51_sym = casos.PS.sym('s51',basis(s51(1)),[length(u) 1]);
s52_sym = casos.PS.sym('s52',basis(s52(1)),[length(u) 1]);

s61_sym = casos.PS.sym('s61',basis(s61(1)),[length(u) 1]);
s62_sym = casos.PS.sym('s62',basis(s62(1)),[length(u) 1]);

s7_sym  = casos.PS.sym('s7',basis(s7));
s8_sym  = casos.PS.sym('s8',basis(s8));
 
sos2 = struct('x',[V;s2;s4;s61;s62;s8], ...
              'p',[V_sym;t0;t1;k;s1_sym;s3_sym;s51_sym;s52_sym;s2_sym;s4_sym;s61_sym;s62_sym]);     

sos2.('g') = [s1_sym*(V) - s2_sym*hT - nabla(V,t) - nabla(V,x)*subs(x_dot,u,k);
              s3_sym.*(V) - s4_sym.*hT - gc;
              subs(u_min,u,k)  + s51_sym.*(V) - s61_sym.*hT;
              subs(u_max,u,k)  + s52_sym.*(V) - s62_sym.*hT;
              s8*(subs(V_sym,t,t0))   - subs(V,t,t0);
              (subs(V,t,t1)) - l];


% states + constraint are SOS cones
opts    = struct;
opts.Kx = struct('s', length(sos2.x(2:end)),'l',1);
opts.Kc = struct('s', length(sos2.g));


S2 = casos.sossol('S2','mosek',sos2,opts);

tend2 = toc;


% profile viewer

%% gamma-beta-V-iteration
disp(['setting up solver 2 took ' num2str(tend2-tend1) ' s'])
disp('==================================================================')

tic
itermax = 10;

t0 = 0;
t1 = T;

disp('Start V-S-iteration ...')
for iter = 1:itermax

    % gamma step
    sol1 = S1('p',[Vval]);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 
 
            disp(['G-step feasible in ' num2str(iter) '/' num2str(itermax)])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['G-step infeasible in ' num2str(iter) '/' num2str(itermax)])
                    return
        % otherwise, error('Failed.')
    end

    % beta step
   sol2 = S2('p',[Vval;t0;t1;sol1.x]);
   switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            Vval = sol2.x(1);

            disp(['B-step feasible in ' num2str(iter) '/' num2str(itermax) ])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['B-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end

end

tend4  = toc;
disp(['Finished G-B-V-iteration after ' num2str(tend4) ])




