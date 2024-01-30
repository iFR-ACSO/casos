% Estimate larget possible safe set


% system states
x = casos.PS('x',2,1);


x_up  = [3 3];
x_low = [-3 -3];

g = [];
for k= 1:length(x_up)
    g = [g ; -(x(k)-x_low(k))*(x_up(k)-x(k) )]; % g(x) <= 0
end

xc = [0;-0];
g =[g ; (x-xc)'*eye(2)*0.1*(x-xc)-1]; % g(x) <= 0


xc = [-1;1];
l = (x-xc)'*eye(2)*2*(x-xc)-1-1e-6; % l(x) > 0

% SOS multiplier
s = casos.PS.sym('s',monomials(x,0:1),[3 1],'gram');
r = casos.PS.sym('r',monomials(x,0:1),'gram');


h = casos.PS.sym('h',monomials(x,0:6));
h_sym = casos.PS.sym('h_sym',basis(h));

% level of stability
c = casos.PS.sym('c');

h_star = x'*eye(2)*10*x-1;

%% Multiplier Step
% profile on
disp('=========================================================')
disp('Build solver...')
tic
% define SOS feasibility
sos = struct('x',[s;r], ...
             'g',[s.*h_sym.*ones(3,1) - g; r*h_sym + l ], ...
             'p',h_sym);

% states + constraint are SOS cones
opts.Kx.s = 4; 
opts.Kx.l = 0; 
opts.Kc.s = 4;
% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S1 = casos.sossol('S1','mosek',sos,opts);

s_sym = casos.PS.sym('s_sym',basis(s(1)),[3 1]);
r_sym = casos.PS.sym('r_sym',basis(r(1)));

% define SOS feasibility
s0 = casos.PS.sym('s0',monomials(x,0:7),'gram');

sos2 = struct('x',[h;s0], ...
              'g',[s_sym.*h.*ones(3,1) - g; r_sym*h + l; s0*h_sym - h ], ...
              'p',[s_sym;r_sym;h_sym]);

% states + constraint are SOS cones
opts.Kx.l = 1; 
opts.Kx.s = 1;
opts.Kc.s = 5;

% ignore infeasibility
opts.error_on_fail = false;

hlb = casos.PS(basis(h),-inf);
hub = casos.PS(basis(h),+inf);

% solve by relaxation to SDP
S2 = casos.sossol('S2','mosek',sos2,opts);
tbuild = toc;
% profile viewer
disp('Finished building solver!')
disp('=========================================================')
disp('Start iteration...')

itermax = 200;
for iter = 1:itermax
    % evaluate parametrized SOS problem
   sol1 = S1('p',h_star);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS'
              % disp(['s step feasible in ' num2str(iter) '/' num2str(itermax) ] )
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
             disp(['s step infeasible in ' num2str(iter) '/' num2str(itermax) ] )
        otherwise, error('Failed.')
    end

     sol2 = S2('p',[sol1.x;h_star],'lbx',hlb,'ubx',hub);

    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS'
             h_star = sol2.x(1);
              % disp(['h step feasible in ' num2str(iter) '/' num2str(itermax) ] )
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
            disp(['h step infeasible in ' num2str(iter) '/' num2str(itermax) ] )
            break
        otherwise, error('Failed.')
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
plotBoxCon([1 2],x_up,x_low)
hold on
pcontour(g(3),0,[-3 3 -3 3],'b')
hold on
pcontour(l,0,[-3 3 -3 3],'r')
pcontour(h_star,0,[-3 3 -3 3],'g')

