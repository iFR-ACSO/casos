% Test and plot the dynamics of a 3-link inverted pendulum

% Number of links
n = 3;
% degree of the polynomial approximation
deg = 2;

% generate 
% 1. full nlink pendulum dynamics
% 2. get the linerize model at the upright position
% 3. get polynomial model of degree 2
[K,S] = generate_nlink_pendulum(n,deg);
save("data.mat", "S","K")

%% Simple Simulation with controlled system
% Simulation time
tspan = [0, 10]; % 10 seconds

% Initial conditions
theta0 = zeros(n,1);   % Initial angles (radians)
theta0(1) = 0.01;
dtheta0 = zeros(n,1);     % Initial angular velocities (rad/s)
initial_state = [theta0; dtheta0];

% Define the ODE function
odefun = @(t, z) pendulum_dynamics_ctrl(z);

% Solve the ODE
[t, z] = ode45(odefun, tspan, initial_state);

% Extract angles for plotting
theta = z(:, 1:n);

% Plot results
figure;
plot(t, theta, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Angle (rad)');
legendStrings = arrayfun(@(idx) sprintf('\\theta_%d', idx), 1:n, 'UniformOutput', false);
legend(legendStrings);
title('Angles of a 3-Link Inverted Pendulum Over Time');
grid on;

%% Get region of attraction

load("data.mat")
% system states
x = casos.PS('x',2*n,1);

% system dynamics
f = pendulum_dynamics_poly(x);

% Lyapunov function candidate
V = x'*S*x;
p = x'*x;

% scalar decision variable
g = casos.PS.sym('g');

% multiplier (idk)
s = casos.PS.sym('s',monomials(x,0:2), 'gram');

% define SOS problem:
sos = struct('x',s,'f',-g,'g', s*(V-g)-nabla(V,x)*f);

% options
opts = struct('sossol','sedumi');

% constraint is scalar SOS cone
opts.Kx = struct('dsos', 1);
opts.Kc = struct('sos', 1);

% solve by relaxation to SDP
S1 = casos.qcsossol('S1','bisection',sos,opts);
%S1 = casos.sossol('S','mosek',sos,opts);
% evaluate
sol = S1();

%%
% import casos.toolboxes.sosopt.*
% import casos.toolboxes.to_multipoly
% figure(1)
% pcontour(subs(V, [x(2:3);x(5:6)], zeros(4,1)), full(sol.x(1)),[-5 5 -5 5])