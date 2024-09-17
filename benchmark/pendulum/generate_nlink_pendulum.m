function [K,S,A_lin,B_lin] = generate_nlink_pendulum(n, deg, only_poly)

% Parameters
m = 1;       % Mass of each link (kg)
l = 1;       % Length of each link (m)
r = 0.05;    % Radius of each link (m)
g = 9.81;    % Gravitational acceleration (m/s^2)

% Moment of inertia for each link (uniform cylindrical rod)
I_center = (1/12) * m * (l^2 + 3*r^2);

% Define symbolic variables
syms theta [1 n] real
syms dtheta [1 n] real
syms ddtheta [1 n] real
syms tau [1 n-1] real

% Position of center of mass for each link
x = sym(zeros(1, n));
y = sym(zeros(1, n));

% Kinetic and Potential Energy
T = 0; % Total Kinetic Energy
V = 0; % Total Potential Energy

for i = 1:n
    % Position of the center of mass of link i
    if i == 1
        x(i) = (l/2)*sin(theta(i));
        y(i) = (l/2)*cos(theta(i));
    else
        x(i) = x(i-1) + (l/2)*sin(theta(i)) + (l/2)*sin(theta(i-1));
        y(i) = y(i-1) + (l/2)*cos(theta(i)) + (l/2)*cos(theta(i-1));
    end
    
    % Velocity of the center of mass of link i
    vx = jacobian(x(i), theta) * dtheta.';
    vy = jacobian(y(i), theta) * dtheta.';
    
    % Kinetic energy (translational + rotational)
    T = T + (1/2)*m*(vx^2 + vy^2) + (1/2)*I_center*dtheta(i)^2;
    
    % Potential energy
    V = V + m*g*y(i);
end

% Lagrangian
L = T - V;

% Equations of Motion (Euler-Lagrange)
eqns = sym(zeros(1, n));
for i = 1:n
    dL_dtheta = diff(L, theta(i));
    dL_ddtheta = diff(L, dtheta(i));
    ddt_dL_ddtheta = jacobian(dL_ddtheta, [theta, dtheta]) * [dtheta, ddtheta].';
    eqns(i) = simplify(ddt_dL_ddtheta - dL_dtheta);
end

% Convert equations to state-space form
f = sym(zeros(2*n, 1)); % State derivatives
z = [theta, dtheta];    % State variables
u = [zeros(1,1), tau];  % Control inputs (tau_1 = 0 for the unactuated joint)

% Resolve system dynamics
% M * ddtheta + C * dtheta + G = B * tau
M = simplify(jacobian(eqns, ddtheta));

AUX = arrayfun(@(idx) diff(M,theta(idx))*dtheta(idx), 1:n, 'UniformOutput', false);
totalSum = 0;
for i=1:n
    totalSum = totalSum + AUX{i};
end

C = jacobian(eqns, dtheta)-0.5*totalSum;
G = subs(eqns, [ddtheta dtheta], zeros(1,2*n));

% Solving for ddtheta
if only_poly == 0
    tic
    ddtheta_sol = vpa(M)*(u' - C*dtheta' - G');
    toc
else
    tic
    % approximate by a polynomial 
    temp_inv = taylor(u' - C*dtheta' - G', [theta, dtheta], 'Order', deg+1, 'ExpansionPoint', zeros(1, 2*n));
    M = taylor(M, [theta, dtheta], 'Order', deg+1, 'ExpansionPoint', zeros(1, 2*n));
    ddtheta_sol = vpa(simplify(M) )*simplify(temp_inv);
    toc
end

% State derivatives
f(1:n) = dtheta';          % dtheta/dt = dtheta
f(n+1:2*n) = ddtheta_sol;  % ddtheta/dt = solution from above
f = (f);

% Convert symbolic to function
pendulum_dyn = ['pendulum_dyn_n' num2str(n)];
matlabFunction(f, 'File', pendulum_dyn, 'Vars', {z', tau'});


%% obtain linear model
% Compute the Jacobian matrices
% Linearized state-space: dx/dt = A*x + B*u
% f is the function handle obtained from the symbolic dynamics
f_sym = feval(pendulum_dyn, [theta, dtheta]', tau');
f_sym = f_sym(:);

% State vector
x_sym = [theta, dtheta];
theta_eq = zeros(1, n);
dtheta_eq = zeros(1, n);
tau_eq = zeros(1, n-1);

% State vector and input vector at equilibrium
x_eq = [theta_eq, dtheta_eq];
u_eq = tau_eq;

% Calculate the Jacobians
A = jacobian(f_sym, x_sym);
B = jacobian(f_sym, tau);

% Substitute equilibrium values into the Jacobian matrices
A_lin = double(subs(A, [x_sym, tau], [x_eq, u_eq]));
B_lin = double(subs(B, [x_sym, tau], [x_eq, u_eq]));

% Define the cost function weights
Q = blkdiag(eye(n), eye(n)); % State cost matrix (weights for the state variables)
R = eye(n-1);        % Control cost matrix (weight for the control inputs)

% Compute the LQR controller gain matrix K
[K, S, e] = lqr(A_lin, B_lin, Q, R);


%% obtain taylor expansion of controlled system
f_sym_ctrl = subs(f_sym, tau', -K*x_sym');
% Convert symbolic to function
matlabFunction(f_sym_ctrl, 'File', ['pendulum_dyn_ctrl_n' num2str(n)], 'Vars', {z'});

taylor_expansion = taylor(f_sym_ctrl, x_sym, 'Order', deg+1, 'ExpansionPoint', x_eq);
% Convert symbolic to function
matlabFunction(taylor_expansion, 'File', ['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)], 'Vars', {z'});

end