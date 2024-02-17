% Estimate larget possible safe set


% system states
x = casos.PS('x',3,1);

% box constraint(s)
x_up  = [ ];
x_low = [];

g = [];

if ~isempty(x_low)
for k= 1:length(x_up)
    g = [g ; -(x(k)-x_low(k))*(x_up(k)-x(k) )]; % g(x) <= 0
end
end

cTheta     = cos(20*pi/180);
epsilon    = 1e-6;

% keep-out
l = [cTheta*(x(1)^2 + x(2)^2 + x(3)^2 + 1)^2 - (8*x(1)*x(2) - x(3)*(4*x(1)^2 + 4*x(2)^2 + 4*x(3)^2 - 4)); ...
     (x(1)^2 + x(2)^2 + x(3)^2 + 1)^2*cTheta + (8*x(1)*x(2) - x(3)*(4*x(1)^2 + 4*x(2)^2 + 4*x(3)^2 - 4))];




% SOS multiplier
s = casos.PS.sym('s',monomials(x,0:1),[length(g) 1],'gram');
r = casos.PS.sym('r',monomials(x,0:1),[length(l) 1],'gram');


% g = [g;l];
s0 = casos.PS.sym('s0',monomials(x,0:2),'gram');

h = casos.PS.sym('h',monomials(x,0:6));
h_sym = casos.PS.sym('h_sym',basis(h));

% level of stability
c = casos.PS.sym('c');

h_star = x'*eye(3)*4*x-0.21;

%% Multiplier Step
% profile on
disp('=========================================================')
disp('Build solver...')
tic
% define SOS feasibility
sos = struct('x',[s;r], ...
             'g',[r.*h_sym.*ones(length(l),1) + l ], ...
             'p',h_sym);

% states + constraint are SOS cones
opts.Kx.s = length(l)+length(g); 
opts.Kx.l = 0; 
opts.Kc.s = length(l)+length(g);

% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S1 = casos.sossol('S1','mosek',sos,opts);

% s_sym = casos.PS.sym('s_sym',basis(s(1)),[length(g) 1]);
r_sym = casos.PS.sym('r_sym',basis(r(1)),[length(l) 1]);

% define SOS feasibility
sos2 = struct('x',[h;s0], ...
              'g',[r_sym.*h.*ones(length(l),1) + l; s0*h_sym - h ], ...
              'p',[r_sym;h_sym]);

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
% profile viewer
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
            break
    end

     sol2 = S2('p',[sol1.x;h_star],'lbx',hlb,'ubx',hub);

    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS'
             h_star = sol2.x(1);
            % figure(1)
            % clf
            %  pcontour3(l(1),0,[-5 5 -5 5 -5 5],'r')
            %     hold on
            %      pcontour3(l(2),0,[-5 5 -5 5 -5 5],'r')
            %  pcontour3(h_star,0,[-5 5 -5 5 -5 5],'r')
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

h_star = cleanpoly(h_star,1e-4);
figure()
for k= 1:length(l)
    pcontour3(l(k),0,[-5 5 -5 5 -5 5],'r')
end
pcontour3(h_star,0,[-5 5 -5 5 -5 5],'g')




