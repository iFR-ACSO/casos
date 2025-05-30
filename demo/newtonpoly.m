% Demonstrate Newton polytope monomial basis reduction 

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

% Enables Newton polytope simplification.
% By default, CaSoS uses the same solver as for the main SDP.
% To disable simplification:
%   - Set opts.newton_solver to []
%   - If using qcsossol, set opts.sossol_options.newton_solver = []
opts.newton_solver = 'mosek'; 

% Build the solver
tic
S = casos.sossol('S','mosek',sos,opts); % build problem
fprintf('build time: %ds \n', toc);

% Run solver
sol = S();
fprintf('sol.f = %d \n', full(sol.f)); % the value should be 1.7107

% Run the solver to obtain an accurate computation time
tsol = timeit(@()S());
fprintf('newton = %s || time: %ds \n', opts.newton_solver, tsol);
