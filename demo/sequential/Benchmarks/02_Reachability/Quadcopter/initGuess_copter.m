
function P = initGuess_copter()

syms x1 x2 x3 x4 x5 x6 u1 u2
x = [x1; x2; x3; x4; x5; x6];
u = [u1; u2];

% parameters
g = 9.81;
K = 0.89/1.4;
d0 = 70;
d1 = 17;
n0 = 55;


% nonlinear dynamics
xdot = [x3; x4; u1*K*sin(x5); ...
        u1*K*cos(x5)-g; x6; ...
        -d0*x5-d1*x6+n0*u2];

Apre = jacobian(xdot,x);
Bpre = jacobian(xdot,u);

% equilibriums
xbar = zeros(6,1);
ubar = [g/K; 0];

% substitute in the value of equilibrium
A = double(subs(subs(Apre,x,xbar),u,ubar));
B = double(subs(subs(Bpre,x,xbar),u,ubar));
Co = ctrb(A,B);

% Q = diag([20 10 20 10 20 10]);
% R = diag([100 100]);
% K = lqr(A,B,Q,R);
% 
% P = lyap((A-B*K)',10*eye(6));

Q = diag([10 10 10 10 10 10]);
R = diag([100 100]);
K = lqr(A,B,Q,R);

P = lyap((A-B*K)',10*eye(6));
