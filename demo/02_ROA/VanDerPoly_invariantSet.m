% Inner-approximation reachable set Van-der-Pol Oscillator

clear
close all
clc
profile off


% indeterminates
x = casos.PS('x',2,1);
u = casos.PS('u',1,1);
w = casos.PS('w');

nx  = length(x);
nu  = length(u);

% system dynamics
f = [x(2);
     (1-x(1)^2)*x(2)-x(1)];

gx = [0;1];

[A,B] = plinearize(f+gx*u,x,u);

[K,P] = lqr(A,B,eye(2),1);

Vval = x'*P*x;

gw = [0;1];

% f = subs(f + gx*u ,u, -K*x);

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials([x;w],2),'gram');
s2 = casos.PS.sym('s2',monomials(x,0)    ,'gram');
s3 = casos.PS.sym('s3',monomials([x;w],2),'gram');
k = casos.PS.sym('s3',monomials(x,0:4));
% enforce positivity
l = 1e-6*(x'*x);

p = x'*x;

% level of stability
g = casos.PS.sym('g');
b = casos.PS.sym('b');

% options
opts               = struct('sossol','mosek');
opts.error_on_fail = 0; 
opts.verbose       = 1;
opts.max_iter      = 20;
opts.conf_interval = [-1000 0]';

wmax = 0.0001;
wcon = (w-(-wmax))*(wmax-w);

%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic

% solver 1: gamma-step
V_sym = casos.PS.sym('V',basis(V));

sos1        = struct('x',[k;s1;s3],'f',-g,'p',V_sym);

sos1.('g')  = s1*(V_sym - g) - nabla(V_sym,x)*(subs(f+gx*u,u,k)+gw*w) - s3*wcon;

% states + constraint 
opts.Kx = struct('s', 2,'l',1);
opts.Kc = struct('s', 1);

S1 = casos.qcsossol('S1','bisection',sos1,opts);

tend1 = toc;
disp(['setting up solver 1 took ' num2str(tend1) ' s'])

% solver 2: beta-step
sos2        = struct('x',s2,'f',-b,'p',[V_sym;g]);
sos2.('g')  = s2*(p-b) + g - V_sym;

% states + constraint are SOS cones
opts.Kx = struct('s', 1);
opts.Kc = struct('s', 1);

S2 = casos.qcsossol('S2','bisection',sos2,opts);

tend2 = toc;
disp(['setting up solver 2 took ' num2str(tend2-tend1) ' s'])

% % solver 3: V-step
% s1_sym = casos.PS.sym('s1',basis(s1));
% s2_sym = casos.PS.sym('s2',basis(s2));
% 
% Vlb = casos.PS(basis(V),-inf);
% Vub = casos.PS(basis(V),+inf);
% 
% sos3       = struct('x',[V;s3],'p',[b,g,s1_sym,s2_sym]);
% sos3.('g') = [V-l; 
%               s2_sym*(p-b)+g-V; 
%               s1_sym*(V-g)-nabla(V,x)*f - s3*(wmax - w'*w)];
% 
% opts    = struct;
% opts.Kx = struct('s', 1, 'l', 1); 
% opts.Kc = struct('s', 3);
% 
% % opts.sdpsol_options.mosek_param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
% % opts.sdpsol_options.mosek_param.MSK_IPAR_INTPNT_BASIS       = 'MSK_BI_NEVER';
% 
% 
% S3 = casos.sossol('S','sedumi',sos3,opts);

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
            s1val = sol1.x(1);
 
            disp(['G-step feasible in ' num2str(iter) '/' num2str(itermax)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['G-step infeasible in ' num2str(iter) '/' num2str(itermax)])
                    
        otherwise, error('Failed.')
    end

    % beta step
    sol2 = S2('p',[Vval;gval]);
   switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 
            bval = double(-sol2.f);
            s2val = sol2.x(1);

            disp(['B-step feasible in ' num2str(iter) '/' num2str(itermax) ])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['B-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end

    % V-step
    sol3 = S3('p',[bval,gval,s1val,s2val],'lbx',Vlb,'ubx',Vub);

   switch (S3.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            Vval = sol3.x(1);
            disp(['V-step feasible in ' num2str(iter) '/' num2str(itermax)])
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end
end

tend4  = toc;
disp(['Finished G-B-V-iteration after ' num2str(tend4) ])

figure()
pcontour(Vval,gval,[-2 2 -2 2])
hold on
pcontour(p,bval,[-2 2 -2 2],'r')
legend('ROA','Shape')












