% Estimate larget possible safe set
clc
clear

% system states
x = casos.PS('x',6,1);

% box constraint(s)


omega_max =  2*pi/180;
x_up  = [ omega_max omega_max omega_max];
x_low = [-omega_max -omega_max -omega_max];


g = [];
for k= 1:length(x_up)
g = [g ;(x(k)-x_low(k))*(x_up(k)-x(k) )]; % g(x) >= 0
end

% keep-in
cTheta = cos(45*pi/180);
epsilon = 1e-6;
g = [g; 
    (x(4)^2 + x(5)^2 + x(6)^2 + 1)^2*(1-cTheta+epsilon) - (8*x(4)^2 + 8*x(6)^2)];

l = [];
% SOS multiplier
s = casos.PS.sym('s',monomials(x,0:1),[length(g) 1],'gram');
% r = casos.PS.sym('r',monomials(x,0:1),[length(l) 1],'gram');

s0 = casos.PS.sym('s0',monomials(x,0:1),'gram');

h = casos.PS.sym('h',monomials(x,0:5));
h_sym = casos.PS.sym('h_sym',basis(h));

% level of stability
% c = casos.PS.sym('c');

h_star = x'*eye(length(x))*2*x;

%% Multiplier Step
profile on
disp('=========================================================')
disp('Build solver...')
tic
% define SOS feasibility
sos = struct('x',s, ...
             'g',s.*h_sym.*ones(length(g),1) + g, ...
             'p',h_sym);

% states + constraint are SOS cones
opts.Kx.s = length(l)+length(g); 
opts.Kx.l = 0; 
opts.Kc.s = length(l)+length(g);

% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S1 = casos.sossol('S1','mosek',sos,opts);

s_sym = casos.PS.sym('s_sym',basis(s(1)),[length(g) 1]);
% r_sym = casos.PS.sym('r_sym',basis(r(1)),[length(l) 1]);

% define SOS feasibility
sos2 = struct('x',[h;s0], ...
              'g',[s_sym.*h.*ones(length(g),1) + g; s0*h_sym - h ], ...
              'p',[s_sym;h_sym]);

% states + constraint are SOS cones
opts.Kx.l = 1; 
opts.Kx.s = 1;
opts.Kc.s = 1+length(g)+length(l);

% ignore infeasibility
opts.error_on_fail = false;

hlb = casos.PS(basis(h),-inf);
hub = casos.PS(basis(h),+inf);

% solve by relaxation to SDP
S2 = casos.sossol('S2','mosek',sos2,opts);
tbuild = toc;


profile viewer
disp('Finished building solver!')
disp('=========================================================')
disp('Start iteration...')

itermax = 100;
for iter = 1:itermax
    % evaluate parametrized SOS problem
   sol1 = S1('p',h_star);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS'
              disp(['s step feasible in ' num2str(iter) '/' num2str(itermax) ] )
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
             disp(['s step infeasible in ' num2str(iter) '/' num2str(itermax) ] )
        otherwise
            disp(['s step infeasible in ' num2str(iter) '/' num2str(itermax) ] )
            % break
    end

     sol2 = S2('p',[sol1.x;h_star],'lbx',hlb,'ubx',hub);

    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS'
             h_star = sol2.x(1);
              disp(['h step feasible in ' num2str(iter) '/' num2str(itermax) ] )
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
            disp(['h step infeasible in ' num2str(iter) '/' num2str(itermax) ] )
            break
        otherwise
            disp(['h step infeasible in ' num2str(iter) '/' num2str(itermax) ] )
            break
    end
end
tIter = toc-tbuild;
disp('=========================================================')
disp('Finished iteration')
disp(['Build time: ' num2str(tbuild) ' s' ])
disp(['Iteration time: ' num2str(tIter) ' s' ])
disp('___________________________________________')
disp(['Total time: ' num2str(tIter+tbuild) ' s' ])

figure()
% plotBoxCon([1 2],x_up,x_low)
% hold on
% 
% for k= 1:length(l)
%     pcontour(l(k),0,[-3 3 -3 3],'r')
% end
% pcontour(h_star,0,[-3 3 -3 3],'g')

% hstruct = to_struct(h_star);
% save('safe_set.mat',"hstruct")

