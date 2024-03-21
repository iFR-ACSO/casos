% Example taken from: 
% Lyapunov Analysis of Rigid Body Systems with Impacts
% and Friction via Sums-of-Squares, Michael Posa, Mark Tobenkin, Russ
% Tedrake

% solve with CaSoS
load hybridLyap.mat;

% decision variables
nX = length(c);
X = casadi.MX.sym('X', nX, 1);

% define SDP problem
sdp.f = c'*X;
sdp.g = A'*X;
sdp.x = X;

% define cones
opts = struct;
opts.Kx = struct('l', K.f ,'s', K.s);

% initialize solver
S = casos.sdpsol('S','sedumi',sdp,opts);

% solve
tic
sol = S('lbg', b, 'ubg', b, 'lbx', -inf, 'ubx', inf);
fprintf('time of computation: %d\n', toc);
