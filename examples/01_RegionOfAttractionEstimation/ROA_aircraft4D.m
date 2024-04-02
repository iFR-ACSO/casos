% ------------------------------------------------------------------------
%
% Short Description: 
% This script implements the the region-of attraction(ROA) estimation 
% problem found in [1]. The resulting bilinear a three step iterative 
% procedure. An elliposoidal shape function is used to ensure grow of
% the sublevel set.
%
% Date: 04/01/2024
%
% Reference: 
% [1] Chakraborty, A., Seiler, P., & Balas, G. J. (2011). 
%  Nonlinear region of attraction analysis for flight control verification 
%  and validation. Control Engineering Practice, 19(4), 335-345. 
%  doi:10.1016/j.conengprac.2010.12.001
%
% ------------------------------------------------------------------------
import casos.toolboxes.sosopt.* % legacy code for plotting


disp(['--------------------------------------------------------------------' ...
    '---------'])
disp(['Compute inner-approximation of the ROA for the longituidnal motion' ...
      ' of the GTM'])
%% Define states and load necessary data
% system states
x = casos.PS('x',4,1);

% Load polynomial dynamics from .mat file; for details see [1]
load GTM_scaled_dyn.mat

f = gtmdyn(x(1),x(2),x(3),x(4));

% shape function
p = x'*x*1e2;

% initial Lyapunov functionn
P = [395.382671347059	-23.0032507976836	3.16965275615691	29.2992065909380
    -23.0032507976836	90.4764915483638	-16.1191428789579	-132.594376986429
    3.16965275615691	-16.1191428789579	3.44002648214490	24.2292666551058
    29.2992065909380	-132.594376986429	24.2292666551058	202.114797577027];

Vval = x'*P*x;

%% Define all polynmoials

% Lyapunov function candidate
V = casos.PS.sym('v',monomials(x,2));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,0:2),'gram');
s2 = casos.PS.sym('s2',monomials(x,0:2),'gram');

% enforce positivity
l = 1e-6*(x'*x);

% level of stability
g = casos.PS.sym('g');
b = casos.PS.sym('b');

% options
opts = struct('sossol','sedumi');

%% setup solvers
disp(['--------------------------------------------------------------------' ...
    '---------'])
disp('Setup solver ...')
tic
% solver 1: gamma-step
prob_gamma_step = struct('x',s1, ...   % decision variables
                         'f',-g, ...   % cost function
                         'p',V);       % parameter

prob_gamma_step.('g') = s1*(V-g)-nabla(V,x)*f-l; % constraint functions

% states + constraint are SOS cones
opts.Kx = struct('sos', 1); % cone for decision variables
opts.Kc = struct('sos', 1); % cone for constraints

% setup solver for quasi-convex problem (biscetion
solver_gamma_step = casos.qcsossol('S1', ...
                                  'bisection', ...
                                   prob_gamma_step, ...
                                   opts);

disp(['... gamma step succesful after ',  num2str(toc) ' s'])
% solver 2: beta-step
tic
prob_beta_step = struct('x',s2, ... % decision variables
                        'f',-b, ... % cost function
                        'p',[V;g]); % parameter

prob_beta_step.('g') = s2*(p-b)+g-V; % constraint functions

% states + constraint are SOS cones
opts.Kx = struct('sos', 1); % cone for decision variables
opts.Kc = struct('sos', 1); % cone for constraints

solver_beta_step = casos.qcsossol('S2', ...
                                  'bisection', ...
                                  prob_beta_step, ...
                                  opts);
disp(['... beta step succesful after ',  num2str(toc) ' s'])
% solver 3: V-step
tic
% parameter; currently only gram-like polynomials can be used for paramter
s1_sym = casos.PS.sym('s1',basis(s1));
s2_sym = casos.PS.sym('s2',basis(s2));

% define lower and upper boubd for Lyapunov function
Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

prob_V_step = struct('x',V, ...                % decision variables
                     'p',[b,g,s1_sym,s2_sym]); % parameter

prob_V_step.('g') = [V-l; 
                     s2_sym*(p-b)+g-V; 
                     s1_sym*(V-g)-nabla(V,x)*f-l]; % constraint functions

% options
opts = struct;
opts.Kx = struct('sos', 0, 'lin', 1);  % cone for decision variables
opts.Kc = struct('sos', 3);            % cone for constraints

solver_V_step = casos.sossol('S', ...
                            'sedumi', ...
                            prob_V_step, ...
                            opts);
disp(['... V step succesful after ',  num2str(toc) ' s'])

%% gamma-beta-V-iteration
disp(['--------------------------------------------------------------------' ...
    '---------'])
disp('Start iterative scheme')
tic
itermax = 10;
for iter = 1:itermax

    % gamma step
    sol_gamma_step = solver_gamma_step('p',Vval); % pass parameter

    switch (solver_gamma_step.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            gval = double(-sol_gamma_step.f);
            s1val = sol_gamma_step.x;
            disp(['G-step feasible in ' num2str(iter) '/' num2str(itermax)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['G-step infeasible in ' num2str(iter) '/' num2str(itermax)])
                    
        otherwise, error('Failed.')
    end

    % beta step
    sol_beta_step = solver_beta_step('p',[Vval;gval]); % pass parameter


   switch (solver_beta_step.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 
            bval = -sol_beta_step.f;
            s2val = sol_beta_step.x;

            disp(['B-step feasible in ' num2str(iter) '/' num2str(itermax) ])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['B-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end


    % V-step
    sol_V_step = solver_V_step('p',[bval,gval,s1val,s2val], ... % pass parameter
                               'lbx',Vlb, ... % set lower bound on dec. variables
                               'ubx',Vub);    % set upper bound on dec. variables

    switch (solver_V_step.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

                Vval = sol_V_step.x;

            disp(['V-step feasible in ' num2str(iter) '/' num2str(itermax) ])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end

end

disp(['End of iterative scheme after: ', num2str(toc),' s'])
disp(['--------------------------------------------------------------------' ...
    '---------'])

%% Plotting
d = ([
          convvel(50, 'm/s', 'm/s') 
          convang(20, 'deg', 'rad')  
          convang(30, 'deg', 'rad')  
          convang(10, 'deg', 'rad') 
]);

D = diag(d)^-1;
V = subs(Vval,x,D*x);
V1 = subs(V,[x(1);x(4)],[0 0]');

p = subs(p,x,D*x);
p1 = subs(p,[x(1);x(4)],[0 0]');

figure(1)
pcontour(V1, double(gval), [-1 1 -4 4], 'b-');
hold on
pcontour(p1, double(bval), [-1 1 -4 4], 'r--');

