clc
clear

pr = spotsosprog;

n = 2; 

[pr,x] = pr.newIndeterminate('x',n);

f = -[ x(2) ; (1-x(1)^2)*x(2)-x(1) ];

A = double(subs(diff(f,x),x,0*x));
P = lyap(A',eye(n));

V = x'*P*x;
Vdot = diff(V,x)*f;

[pr,r] = pr.newFree(1);
[pr,l] = pr.newFreePoly(monomials(x,0:4));

d = floor(deg(l*Vdot,x)/2-1);
pr = pr.withSDSOS(l*Vdot + (V-r)*(x'*x)^d);

opt = spot_sdp_default_options();
opt.verbose = 0;

tic
sol = pr.minimize(-r,@spot_mosek,opt);
toc

% Optimal value
opt_sos = double(sol.eval(r));

disp(['Optimal value (SDSOS): ' num2str(opt_sos)])


pr = pr.withSOS(l*Vdot + (V-r)*(x'*x)^d);

opt = spot_sdp_default_options();
opt.verbose = 0;

tic
sol = pr.minimize(-r,@spot_mosek,opt);
toc

% Optimal value
opt_sos = double(sol.eval(r));

disp(['Optimal value (SOS): ' num2str(opt_sos)])