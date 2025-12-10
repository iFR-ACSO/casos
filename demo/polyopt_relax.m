% Polynomial optimization.
% Use of smaller cones such as DSOS and SDSOS
%% Sum-of-squares cone
% indeterminate variable
x = casos.Indeterminates('x');
% some polynomial
f = x^4 + 1*x;
% scalar decision variable
g = casos.PS.sym('g');

% define SOS problem:
%   min g s.t. (f + g) is SOS
sos = struct('x',g,'f',g,'g',f+g);

%% Sum-of-squares cone

% constraint is scalar SOS cone
opts = struct('Kc',struct('sos',1));

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);
% evaluate
tic
sol = S();
% result should be 0.47247
fprintf('SOS: Minimum is %g. (time = %d)\n', full(sol.f), toc)
% save solution with SOS cone
g_sos = full(sol.x);

%% Diagonally dominant sum-of-squares cone

% constraint is scalar DSOS cone
opts = struct('Kc',struct('dsos',1));

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);
% evaluate
tic
sol = S();
% the result should be 0.75 
fprintf('DSOS: Minimum is %g. (time = %d)\n', full(sol.f), toc)
% save solution with DSOS cone
g_dsos = full(sol.x);

%% Scaled diagonally dominant sum-of-squares cone

% constraint is scalar SDSOS cone
opts = struct('Kc',struct('sdsos',1)); 

% solve by relaxation to SDP
S = casos.sossol('S','mosek',sos,opts);
% evaluate
tic
sol = S();
% the result should be 0.47247
fprintf('SDSOS: Minimum is %g. (time = %d)\n', full(sol.f), toc)
% save solution with SOS cone
g_sdsos = full(sol.x);

%% Plot results
% convert f into a function to evaluate on the interval [-1, 1]
f_fun = f.to_function;
grid_x = linspace(-1, 1);

% setup figure 
figure(1)
xlabel('x')
ylabel('y')
hold on
legend('show');
grid 
grid minor

% plot the horizontal line with zero
plot(grid_x, zeros(length(grid_x),1), '--', 'DisplayName','Zero line')

% plot f
plot(grid_x, full(f_fun(grid_x)), 'DisplayName','Original polynomial f')

% plot f+g with the g that was found with SOS, DSOS and SDSOS
plot(grid_x, full(f_fun(grid_x)) + g_sos, 'DisplayName','f+g (with SOS)')
plot(grid_x, full(f_fun(grid_x)) + g_dsos, 'DisplayName','f+g (with DSOS)')
plot(grid_x, full(f_fun(grid_x)) + g_sdsos, 'DisplayName','f+g (with SDSOS)')
