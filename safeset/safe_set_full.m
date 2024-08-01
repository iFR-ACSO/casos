% Estimate larget possible safe set
clear
close all
clc


import casos.toolboxes.sosopt.*
% system states
x = casos.PS('x',6,1);

% box constraint(s)
x_up  = [3 3 3 ];
x_low = -x_up;

g = [];

if ~isempty(x_low)
for k= 1:length(x_up)
    g = [g ; -(x(k+3)-x_low(k))*(x_up(k)-x(k+3) )]; % g(x) <= 0
end
end

cTheta     = cos(20*pi/180);
epsilon    = 1e-6;

% keep-out
% l = [cTheta*(x(1)^2 + x(2)^2 + x(3)^2 + 1)^2 - (8*x(1)*x(2) - x(3)*(4*x(1)^2 + 4*x(2)^2 + 4*x(3)^2 - 4)); ...
%      cTheta*(x(1)^2 + x(2)^2 + x(3)^2 + 1)^2 + (8*x(1)*x(2) - x(3)*(4*x(1)^2 + 4*x(2)^2 + 4*x(3)^2 - 4))];


% l1 < 0  --> fulfilled

% l2 < 0 --> fullfilled

% l1 / l2 >= 0
l = [cTheta*(x(1)^2 + x(2)^2 + x(3)^2 + 1)^2 - (8*x(1)*x(2) - x(3)*(4*x(1)^2 + 4*x(2)^2 + 4*x(3)^2 - 4));  
      (x(1)^2 + x(2)^2 + x(3)^2 + 1)^2*cTheta + (8*x(1)*x(2) - x(3)*(4*x(1)^2 + 4*x(2)^2 + 4*x(3)^2 - 4))];

% SOS multiplier
s = casos.PS.sym('s',monomials(x,0:1),[length(g) 1],'gram');
r = casos.PS.sym('r',monomials(x,0:1),[length(l) 1],'gram');

% g = [g;l];
s0 = casos.PS.sym('s0',monomials(x,0:2),'gram');

h = casos.PS.sym('h',monomials(x,0:8));
h_sym = casos.PS.sym('h_sym',sparsity(h));

% level of stability
c = casos.PS.sym('c');

h_star = x'*eye(6)*5*x-0.01;


figure()
pcontour3(subs(h_star,x(4:6),zeros(3,1)),0,[-1 1 -1 1 -1 1])

%% Multiplier Step
% profile on
disp('=========================================================')
disp('Build solver...')
tic
% define SOS feasibility
sos = struct('x',[s;r], ...
             'g',[r*h_sym + l - epsilon; ...
                  s*h_sym - g], ...
             'p',h_sym);

% states + constraint are SOS cones
opts.Kx.sos = length(l)+length(g); 
opts.Kx.lin = 0; 
opts.Kc.sos = length(l)+length(g);

% ignore infeasibility
opts.error_on_fail = false;

% solve by relaxation to SDP
S1 = casos.sossol('S1','mosek',sos,opts);

s_sym = casos.PS.sym('s_sym',sparsity(s(1)),[length(g) 1]);
r_sym = casos.PS.sym('r_sym',sparsity(r(1)),[length(l) 1]);

% define SOS feasibility
sos2 = struct('x',[h;s0], ...
              'g',[r_sym*h + l - epsilon; ...
                   s_sym*h - g;
                   s0*h_sym - h ], ...
              'p',[r_sym;s_sym;h_sym]);

% states + constraint are SOS cones
opts.Kx.lin = 1; 
opts.Kx.sos = 1;
opts.Kc.sos = 1+length(g)+length(l);

% ignore infeasibility
opts.error_on_fail = false;

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

     sol2 = S2('p',[sol1.x;h_star]);

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

h_star = cleanpoly(h_star,1e-4);

%% keep out constraints
import casos.toolboxes.sosopt.*
figure()
for k= 1:length(l)
    pcontour3(1*casos.PD(l(k)),0,[-5 5 -5 5 -5 5])
end

% plot safe set as 3D volume
pcontour3(subs(h_star,x(4:6),zeros(3,1)),0,[-5 5 -5 5 -5 5],'g')

%% plot safe set on unit sphere

% Define direction vector to star w.r.t. to the inertial frame 
starVec1_I = [0;1;0];
starVec2_I = [0;-1;0];

cone_angle1 = 20*pi/180;
cone_angle2 = 20*pi/180;

% Define boresight vector of instrument in body frame;
borSight_B = [1;0;0];

[x,y,z] = sphere(100);

figure()
mesh(x,y,z,'EdgeColor',[.4 .4 .4],'FaceAlpha',0.0,'EdgeAlpha',0.1)
axis equal
hold on
% plot coordinate system
plot3([0 0 1], [0 0 0], [0 0 0], 'r','LineWidth',2)
plot3([0 0 0], [0 0 1], [0 0 0], 'g','LineWidth',2)
plot3([0 0 0], [0 0 0], [0 0 1], 'b','LineWidth',2)
% plot cones
PlotCone(starVec1_I, cone_angle1,[1 0 0]);
PlotCone(starVec2_I, cone_angle2,[1 0 0]);
box on

% generate attitude in MRP space
a = -1;
b = 1;

sigma = a +(b-a)*rand(3,10000);
sigma = sigma./vecnorm(sigma);


a = -3*pi/180;
b = 3*pi/180;

omega = a +(b-a)*rand(3,10000);

% get only attitude that lie in safe set
hfun         = to_function(h_star);
hfun_values  = full(hfun(sigma(1,:),sigma(2,:),sigma(3,:),omega(1,:),omega(2,:),omega(3,:)));


idx = find(hfun_values <= 0);

% plot safe attitudes in unit-sphere space i.e. we want that the borsight
% angle stays outside of the keep-out cones
for k = 1:length(idx)

    T_IB = mrp2trafo(sigma(:,idx(k)))';
    borSight_I = T_IB*borSight_B;
    plot3(borSight_I(1),borSight_I(2),borSight_I(3),'*g')

end





