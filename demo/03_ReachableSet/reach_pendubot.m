% ------------------------------------------------------------------------
%
%   Short Description:  Inner-approximation of the reachable set of the GTM
%                       model with respect to control constraint and
%                       terminal set.
%
%   Source: Yin et al. 
%
%   Date: 02/17/2024
%
% ------------------------------------------------------------------------


clear
close all
clc


%% define variables
x = casos.PS('x',4,1);
u = casos.PS('u',1,1);

% time
t  = casos.PS('t');

x1 = x(1);
x2 = x(2); 
x3 = x(3); 
x4 = x(4); 


%% define your system here

% scaled 4-state system
f2 = - 10.6560*x1^3 + 11.5309*x1^2*x3 + 7.8850*x1*x3^2 + 0.7972*x2^2*x3 ...
  + 0.8408*x2*x3*x4 + 21.0492*x3^3 + 0.4204*x3*x4^2 + 66.5225*x1 - 24.5110*x3;

f4 = 10.9955*x1^3 - 48.9151*x1^2*x3 - 6.4044*x1*x3^2 - 2.3955*x2^2*x3 ...
  - 1.5943*x2*x3*x4 - 51.9088*x3^3 - 0.7971*x3*x4^2 - 68.6419*x1 + 103.9783*x3;

f = [x2; f2; x4; f4];

g2 = -10.0959*x3^2 + 44.2521;
g4 =  37.8015*x3^2 - 83.9120;

gx = [0; g2; 0; g4];

% dynamics
x_dot =  f+gx*u;

% target function
l = x'*blkdiag(1/0.1^2, 1/0.35^2, 1/0.1^2, 1/0.35^2)*x - 1;

% l =x'*x-1;
Vval = 182.2102*x1^2 + 67.8981*x1*x2 + 314.3265*x1*x3 + 37.2705*x1*x4 + ...
    6.4123*x2^2 + 59.0528*x2*x3 + 7.0314*x2*x4 + 138.9343*x3^2 + ...
    32.5944*x3*x4 + 1.9521*x4^2;


figure()
pcontour(subs(l,x(3:4),zeros(2,1)),0,[-1 1 -1 1])
hold on
pcontour(subs(Vval,x(3:4),zeros(2,1)),0.1,[-1 1 -1 1])

figure()
pcontour(subs(l,x(1:2),zeros(2,1)),0,[-1 1 -1 1])
hold on
pcontour(subs(Vval,x(1:2),zeros(2,1)),0.1,[-1 1 -1 1])

% Trim point for elevator channel is 0.0489 rad.
% Saturation limit for elevator channel is -10 deg to 10 deg
um =  -1; 
uM =   1;

% affine control constraints
u_min = u - um; 
u_max = uM - u;


%% define SOS problem polynomials

% reachability storage function
V = casos.PS.sym('v',monomials([x;t],0:2));
k = casos.PS.sym('k',monomials([x;t],0:2));

% SOS multiplier dissipation inequality
s1 = casos.PS.sym('s1',monomials([x;t],0:2));
s2 = casos.PS.sym('s2',monomials([x;t],0:2));

% SOS multiplier for control constraints
s51 = casos.PS.sym('s51',monomials([x;t],0:2));
s52 = casos.PS.sym('s52',monomials([x;t],0:2));

s61 = casos.PS.sym('s61',monomials([x;t],0:2));
s62 = casos.PS.sym('s62',monomials([x;t],0:2));

% SOS multiplier terminal-set inclusion
s4 = casos.PS.sym('s4',monomials(x,0:2));  

g = casos.PS.sym('g');

% Time horizon
T  = 4;

t0 = 0;
t1 = T;

% time polynomial
hT = (t)*(T-t);


%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic

% V_sym = casos.PS.sym('V',basis(V));

sos1 = struct('x',[k;s2;s1;s51;s52;s61;s62;s4],...
              'p',[V;g]);     

sos1.('g') = [s2;s1;s51;s52;s61;s62;
              % s4-0.001;
              s1*(V-g) - s2*hT - nabla(V,t) - nabla(V,x)*(f+gx*k);
              k - um     + s51.*(V-g) - s61*hT;
              uM - k     + s52.*(V-g) - s62*hT;
              % subs(V,t,T) - g  - s4*l
              ];

% states + constraint are SOS cones
% opts = struct('sossol','mosek');
opts.Kx = struct('s', size(sos1.x( (length(u)+1):end) ,1), ...
                 'l',length(u));                

opts.Kx = struct('s', 0, ...
                 'l',size(sos1.x,1));                


opts.Kc = struct('s', length(sos1.g));          % we only have sos con.

% setupsolve
opts.error_on_fail = 0; 
% opts.verbose = 1;
% opts.max_iter      = 25;
% opts.conf_interval = [-1 0]';
% opts.tolerance_abs = 1e-4;
% opts.tolerance_rel = 1e-4;

solve_Sstep = casos.sossol('S1','mosek',sos1,opts);

tend1 = toc;

disp(['setting up solver 1 took ' num2str(tend1) ' s'])



%% V-s-iteration
% disp(['setting up solver 2 took ' num2str(tend2) ' s'])

disp('==================================================================')

tic

itermax = 10;

disp('Start V-S-iteration ...')
for iter = 1:itermax
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
        gamma_try = 0.015625000000000;
        % find mulitplier and control law
        sol_Sstep = solve_Sstep('p',[Vval;gamma_try]);
       solve_Sstep.stats.UNIFIED_RETURN_STATUS
        switch (solve_Sstep.stats.UNIFIED_RETURN_STATUS)
            case 'SOLVER_RET_SUCCESS' 
                 
                 gamma_lb = gamma_try;
                 go = 0;
                disp(['S-step feasible in ' num2str(iter) '/' num2str(itermax)])
                
            case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
                gamma_ub = gamma_try;
                disp(['g-step infeasible in ' num2str(iter) '/' num2str(itermax)])
          
        end
    end

   % find reachability storage function
   sol_Vstep = solver_Vstep('p',[-sol_Sstep.f;Vval;sol_Sstep.x(1:7)]);


   switch (solver_Vstep.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            Vval = sol_Vstep.x(1);

            disp(['V-step feasible in ' num2str(iter) '/' num2str(itermax) ])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end

end

tend4  = toc;
disp(['Finished V-s-iteration after ' num2str(tend4) ])




