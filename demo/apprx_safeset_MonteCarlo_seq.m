% Estimate largest possible safe set
clear 
close all
clc

import casos.toolboxes.sosopt.*

% system states
x = casos.PS('x',2,1);

% box constraint(s)
x_up  = [ 3 ; 3];
x_low = [-3; -3];

g = -(x-x_low).*(x_up-x);

% keep-out
xc = [-2;2];
l = (x-xc)'*eye(2)*2*(x-xc)-1-1e-6; % l(x) > 0

xc = [2;-2];
l =[l; (x-xc)'*eye(2)*2*(x-xc)-1-1e-6]; % l(x) > 0

% l = [];

% SOS multiplier
s = casos.PS.sym('s',monomials(x,0),[length(g) 1]);
r = casos.PS.sym('r',monomials(x,0),[length(l) 1]);
b = casos.PS.sym('b');
h = casos.PS.sym('h',monomials(x,0:6));

% superellipsoid
g0 = 0;
n0 = 50;
for k = 1:2

    g0  = g0 + (x(k)^2/x_up(k)^2)^(n0/2) ;

end
g0 = g0-1;

% plot problem
figure(1)
pcontour(g0,0,[-4 4 -4 4],'k')
hold on
pcontour(l(1),0,[-3 3 -3 3],'r')
pcontour(l(2),0,[-3 3 -3 3],'r')


disp('=========================================================')
disp('Build solver...')
tic

opts = struct('sossol','mosek');


cost_fun = to_function( h );
% keep_out_fun = to_function(l);
% keep_in_fun = to_function(g0);
% 
lower = -2.3;
upper =  2.3;
% 
p = (upper-lower)*net(haltonset(2),100)'+lower;
% 
% l1_val = full(casadi.DM(keep_out_fun(casadi.SX(p))));
% g1_val = full(casadi.DM(keep_in_fun(casadi.SX(p))));
% 
% idx = find(all(([( logical(l1_val > 0) )'  (logical(g1_val <= 0) )']),2));
% 
% 
% plot(p(1,idx),p(2,idx),'k*')




cost = sum(cost_fun(casadi.SX(p)))/100;

% problem struct
sos = struct('x',[h;s;r;b], ...
             'f',cost,...
             'g',[s; r; s*(h-b) - g; r*(h-b) + l  ]);

% states + constraint are SOS cones
opts.Kx.sos = 0; 
opts.Kx.lin = length(sos.x); 
opts.Kc.sos = length(sos.g);
opts.Sequential_Algorithm = 'SQP';
opts.verbose = 1;

% upper and lower bounds for decision variables i.e. coefficients
hlb = casos.PS(basis(h), -inf);
hub = casos.PS(basis(h), +inf);
slb = casos.PS(basis(s), -inf);
sub = casos.PS(basis(s), +inf);
rlb = casos.PS(basis(r), -inf);
rub = casos.PS(basis(r), +inf);
blb = casos.PS(basis(b), -inf);
bub = casos.PS(basis(b), +inf);

lbx = [hlb;slb;rlb;blb];
ubx = [hub;sub;rub;bub];

% setup solver
tic
S1 = casos.nlsossol('S1','sequential',sos,opts);
toc

h0 = x'*eye(2)*x-1;
s0 = [1;1];
r0 = [1;1];
b0 = 1;

% solve problem
sol = S1('x0' ,[h0;s0;r0;b0], ...
         'lbx',lbx, ...
         'ubx',ubx);

% plot beta level set
pcontour(sol.x(1),double(sol.x(end)),[-3 3 -3 3],'g--')
