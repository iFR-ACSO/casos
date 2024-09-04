% Test Newton polytope monomial basis reduction 

x = casos.PS('x', 3, 1);
f = casos.PS.sym('f');

g = [12+x(2)^2-2*x(1)^3*x(2)+2*x(2)*x(3)^2+x(1)^6-2*x(1)^3*x(3)^2+x(3)^4+x(1)^2*x(2)^2+f;
    (x(1)-1)^2+x(1)^4+f-2;
    3*x(1)^4-2*x(1)^2*x(2)+7*x(1)^2-4*x(1)*x(2)+4*x(2)^2+f;
    x(1)^4 + x(2)^4 - 2*x(1)^2*x(2)^2];

sos = struct('x', f, 'g', g, 'f', f);
% states + constraint are SOS cones
opts.Kx.lin = 1;
opts.Kx.sos = 0;
opts.Kc.sos = 4;

% ignore infeasibility
opts.error_on_fail = false;

% newton polytope simplification 
opts.sossol_options.newton = 1;

% solve by relaxation to SDP
tic
S = casos.sossol('S','mosek',sos,opts); % build problem
fprintf('build time: %ds \n', toc);

% to easily verify the computation time, set number of repetitions
N_times = 1000;
idx = 1:N_times;

tic;
arrayfun(@(i) S(), idx, 'UniformOutput', false); % solve problem
fprintf('newton = %d || time: %ds \n', opts.sossol_options.newton, toc);

sol = S();
fprintf('sol.f = %d \n', full(sol.f)); % the value should be 1.7107





