% Test and plot the dynamics of a 3-link inverted pendulum

% Number of links
n = 2;

% degree of the polynomial approximation
deg = 2;

pendulum_dyn_ctrl = ['pendulum_dyn_ctrl_n' num2str(n)];
pendulum_dyn_poly = ['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)];

% generate 
% 1. full nlink pendulum dynamics
% 2. get the linerize model at the upright position
% 3. get polynomial model of degree 2
data_filename = ['data_n' int2str(n) '.mat'];
if exist(data_filename, 'file')~=2 || ...
        exist(pendulum_dyn_poly, 'file')~=2
    [K,S] = generate_nlink_pendulum(n,deg);
    save(data_filename, "S","K")
else
    load(data_filename)
end

%% Get region of attraction

% system states
x = casos.PS('x',2*n,1);

% system dynamics
f = feval(['pendulum_dyn_poly_n' num2str(n) '_d' num2str(deg)],x);

% Lyapunov function candidate
V = x'*S*x;
p = x'*x;

% scalar decision variable
g = casos.PS.sym('g');

% multiplier (idk)
s = casos.PS.sym('s',monomials(x,0:2));

% define SOS problem:
sos = struct('x',s,'f',-g,'g', s*(V-g)-nabla(V,x)*f);

% options
opts = struct('sossol','mosek');

% constraint is scalar SOS cone
opts.Kx = struct('sos', 1);
opts.Kc = struct('sos', 1);

% solve by relaxation to SDP
S1 = casos.qcsossol('S1','bisection',sos,opts);

% evaluate
tic
sol = S1();
fprintf('%s:%d\n', S1.stats.UNIFIED_RETURN_STATUS, toc);
%disp(S1.stats.UNIFIED_RETURN_STATUS)

%%
import casos.toolboxes.sosopt.*
import casos.toolboxes.to_multipoly
figure(1)
V_2d = subs(V, [x(2:n); x(n+2:2*n)], zeros(2*n-2,1));
level = full(-sol.f(1));
pcontour(V_2d, level, [-5 5 -12 12])
grid 
grid minor


