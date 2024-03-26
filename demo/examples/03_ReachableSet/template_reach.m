clc
clear
close all
profile off

%% define variables
x = casos.PS('x',[],1);
u = casos.PS('u',[],1);

% time
t  = casos.PS('t');
t0 = casos.PS.sym('t0');
t1 = casos.PS.sym('t1');



%% define your system here

% dynamics
x_dot =  []


% gc expects functions of the form g(x) <= 0
gc = []

umin = [] 
umax = []

% affine control constraints
u_min = u - umin; 
u_max = umax - u;


%% define SOS problem polynomials

% reachability storage function
V = casos.PS.sym('v',monomials([x;t],));
k = casos.PS.sym('k',monomials([x;t],),[length(u) 1]);

% SOS multiplier dissipation inequality
s1 = casos.PS.sym('s1',monomials([x;t],),'gram');
s2 = casos.PS.sym('s2',monomials([x;t],),'gram');

% SOS multiplier for state constraints
s3 = casos.PS.sym('s3',monomials([x;t],),[length(gc) 1],'gram');
s4 = casos.PS.sym('s4',monomials([x;t],),[length(gc) 1],'gram');

% SOS multiplier for control constraints
s51 = casos.PS.sym('s51',monomials([x;t],),[length(u) 1],'gram');
s52 = casos.PS.sym('s52',monomials([x;t],),[length(u) 1],'gram');

s61 = casos.PS.sym('s61',monomials([x;t],),[length(u) 1],'gram');
s62 = casos.PS.sym('s62',monomials([x;t],),[length(u) 1],'gram');

% SOS multiplier terminal-set inclusion
s7 = casos.PS.sym('s7',monomials(x,0),'gram');  

% SOS multiplier grow constraint
s8 = casos.PS.sym('s8',monomials(x,0),'gram');   

% Time horizon
T  = []

t0 = 0;
t1 = T;

% time polynomial
hT = (t-t0)*(t1-t);


%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic

V_sym = casos.PS.sym('V',basis(V));

sos1 = struct('x',[k;s1;s3;s51;s52;s2;s4;s61;s62],...
              'p',V_sym);     

sos1.('g') = [s1*(V_sym) - s2*hT - nabla((V_sym),t) - nabla((V_sym),x)*subs(x_dot,u,k);
              s3.*(V_sym) - s4.*hT  - gc;
              subs(u_min,u,k)     + s51.*(V_sym) - s61.*hT;
              subs(u_max,u,k)     + s52.*(V_sym) - s62.*hT];

% states + constraint are SOS cones
opts    = struct;
opts.Kx = struct('s', size(sos1.x( (length(u)+1):end) ,1), ...
                 'l',length(u));                

opts.Kc = struct('s', length(sos1.g));          % we only have sos con.

% setupsolve
solve_Sstep = casos.sossol('S1','mosek',sos1,opts);

tend1 = toc;

disp(['setting up solver 1 took ' num2str(tend1) ' s'])

% copy SOS constraints; currently SOS poly cannot be used for
% parameterization; and helps code readability
s1_sym = casos.PS.sym('s1',basis(s1));
s2_sym = casos.PS.sym('s2',basis(s1));
s3_sym = casos.PS.sym('s3',basis(s3(1)),[length(gc) 1]);
s4_sym = casos.PS.sym('s4',basis(s4(1)),[length(gc) 1]);

s51_sym = casos.PS.sym('s51',basis(s51(1)),[length(u) 1]);
s52_sym = casos.PS.sym('s52',basis(s52(1)),[length(u) 1]);

s61_sym = casos.PS.sym('s61',basis(s61(1)),[length(u) 1]);
s62_sym = casos.PS.sym('s62',basis(s62(1)),[length(u) 1]);

s7_sym  = casos.PS.sym('s7',basis(s7));
s8_sym  = casos.PS.sym('s8',basis(s8));

k_sym = casos.PS.sym('V',basis(k));

sos2 = struct('x',[V;s8], ...
              'p',[V_sym;k_sym;s1_sym;s3_sym;s51_sym;s52_sym;s2_sym;s4_sym;s61_sym;s62_sym]);     

sos2.('g') = [
              s1_sym*(V) - s2_sym*hT - nabla(V,t) - nabla(V,x)*subs(x_dot,u,k_sym);
              subs(u_min,u,k_sym)  + s51_sym.*(V) - s61_sym.*hT;
              subs(u_max,u,k_sym)  + s52_sym.*(V) - s62_sym.*hT;
              s3_sym.*V - s4_sym.*hT - gc;
              subs(V,t,t1) - l
              s8*(subs(V_sym,t,t0))   - subs(V,t,t0);
              ];


% states + constraint are SOS cones
opts    = struct;
opts.Kx = struct('s', length(sos2.x(2:end)), ... % first entry storage fun.
                 'l',1); 
opts.Kc = struct('s', length(sos2.g));           % we only have sos con.

% setup solver
solver_Vstep = casos.sossol('S2','mosek',sos2,opts);

tend2 = toc;



%% V-s-iteration
disp(['setting up solver 2 took ' num2str(tend2-tend1) ' s'])

disp('==================================================================')

tic

itermax = 10;

disp('Start V-S-iteration ...')
for iter = 1:itermax

    % find mulitplier and control law
    sol_Sstep = solve_Sstep('p',Vval);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 
 
            disp(['S-step feasible in ' num2str(iter) '/' num2str(itermax)])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter) '/' num2str(itermax)])
                    return

    end

   % find reachability storage function
   sol_Vstep = solver_Vstep('p',[Vval;sol_Sstep.x]);


   switch (S2.stats.UNIFIED_RETURN_STATUS)
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




