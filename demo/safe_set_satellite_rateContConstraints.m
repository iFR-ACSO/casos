% Estimate largest possible safe set
clear 
close all
clc

import casos.toolboxes.sosopt.*

% system states
x = casos.PS('x',3,1);
u = casos.PS('u',3,1);

omega_max = 2*pi/180;

% box constraint(s)
x_up  = [ omega_max ; omega_max; omega_max];
x_low = [-omega_max; -omega_max; -omega_max];

g1 = -(x-x_low).*(x_up-x);


% Casini Parameter
J = diag([8970;9230;3830]);

% skew-symmetric matrix
skew = @(x) [   0  -x(3)  x(2); 
              x(3)   0   -x(1); 
             -x(2)  x(1)   0 ];

% system dynamics
f = -inv(J)*skew(x(1:3))*J*x(1:3);

gx = inv(J);

[A,B] = plinearize(f + gx*u ,x , u);

[K0,~] = lqr(A,B,eye(3)*0.1,eye(3)*0.01);

K = -K0*x;

u_up  =  ones(3,1).*0.1; % Nm or 100 mNm
u_low = -ones(3,1).*0.1;


% superellipsoid
g0 = 0;
n0 = 8;
for k = 1:3

    g0  = g0 + (x(k)^2/x_up(k)^2)^(n0/2) ;

end
g0 = g0-1;



g2 = 0;
for k = 1:3

    g2  = g2 + (K(k)^2/u_up(k)^2)^(n0/2) ;

end
g2 = g2-1;


g = cleanpoly([g0;g2],1e-5);

l = [];

% SOS multiplier
s = casos.PS.sym('s',monomials(x,0),[length(g) 1]);
r = casos.PS.sym('r',monomials(x,0),[length(l) 1]);
b = casos.PS.sym('b');
h = casos.PS.sym('h',monomials(x,0:n0));

disp('=========================================================')
disp('Build solver...')
tic

opts = struct('sossol','mosek');


cost = dot(g0 - (h-b), g0 - (h-b));

% problem struct
sos = struct('x',[h;s;r;b], ...
             'f',cost,...
             'g',[s; r; s*(h-b) - g; r*(h-b) + l  ]);

% states + constraint are SOS cones
opts.Kx.sos = 0; 
opts.Kx.lin = length(sos.x); 
opts.Kc.sos = length(sos.g);
opts.Sequential_Algorithm = 'SQP';
opts.verbose                = 1;

% upper and lower bounds for decision variables i.e. coefficients
hlb = casos.PS(basis(h), -inf);
hub = casos.PS(basis(h), +inf);
slb = casos.PS(basis(s), -inf);
sub = casos.PS(basis(s), +inf);
rlb = casos.PS(basis(r), -inf);
rub = casos.PS(basis(r), +inf);
blb = casos.PS(basis(b), 0);
bub = casos.PS(basis(b), +inf);

lbx = [hlb;slb;rlb;blb];
ubx = [hub;sub;rub;bub];

% setup solver
tic
S1 = casos.nlsossol('S1','sequential',sos,opts);
toc

h0 = g0;
s0 = [1;1];
r0 = [];
b0 = 0;


figure(1)
pcontour(subs(h0,x(3),0),b0,[-omega_max omega_max -omega_max omega_max],'b')
hold on
pcontour(subs(g0,x(3),0),0,[-omega_max omega_max -omega_max omega_max],'k--')

% solve problem
sol = S1('x0' ,[h0;s0;r0;b0], ...
         'lbx',lbx, ...
         'ubx',ubx);


figure(1)
pcontour(subs(sol.x(1),x(3),0),double(sol.x(end)),[-omega_max omega_max -omega_max omega_max],'g')
hold on
pcontour(subs(g0,x(3),0),0,[-omega_max omega_max -omega_max omega_max],'k--')

h_sol = cleanpoly(sol.x(1)-double(sol.x(end)),1e-10)
