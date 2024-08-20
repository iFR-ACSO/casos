%---------------------------------------------------------------------
% gsosoptdemo1
%
% Demonstration of GSOSOPT function for solving quasisos optimizations.
% This example uses GSOSOPT to compute an estimate of the region of 
% attraction for the van der Pol (VDP) oscillator using the Lyapunov 
% function obtained via linearization. See GSOSOPT help for more details 
% on the function syntax.
%
% Reference: 
% P. Seiler and G. Balas, Quasiconvex Sum-of-Squares Programming,
% IEEE Conference on Decision and Control, 2010.
%---------------------------------------------------------------------

% Create vector field for VDP dynamics: dx/dt = f(x)
pvar x1 x2;
x = [x1;x2];
x1dot = -x2;
x2dot = x1+(x1^2-1)*x2;
f = [x1dot; x2dot];

% Plot VDP limit cycle using backward simulation. The actual region of
% attraction is the interior of the set bounded by this limit cycle.
X0 = [1/2; 1/2];
Tf = 100;
[xtraj,isconv]=psim(-f,x,X0,Tf);
t = xtraj{1};
xtraj = xtraj{2};
idx = min(find(t>=Tf*0.8));

figure(1)
hold off;
ph1=plot(xtraj(idx:end,1),xtraj(idx:end,2));
set(ph1,'LineWidth',2)
xlabel('x1');
ylabel('x2');
hold on;

% Construct Lyapunov function from linearization 
Alin = plinearize(f,x);
Q = eye(2);
P=lyap(Alin',Q);
V = x'*P*x;

% Define monomials for S-procedure multiplier, s(x)
z = monomials(x, 1:2 );

% Estimate Region of Attraction (ROA)
% If { x : V(x)<=gamma } \in { x : dV/dt < 0 OR x=0 } then
% { x : V(x)<=gamma } is a subset of the region of attraction.
% To find the largest estimate of the region of attraction using V(x),
% solve the quasisos problem:
%
% max gamma 
% subject to:
% { x : V(x)<=gamma } \in { x : dV/dt < 0 OR x=0 }
%
% A bilinear SOS sufficient condition can be formulated for the set 
% containment constraint and then solved using gsosopt.

% Set options
opts = gsosoptions;
opts.minobj = -50; 
opts.maxobj = 0;


% Form bilinear constraints.  The variable t:=-gamma is used to convert
% the maximization of gamma into a minimization of t.
pvar t;
s = sosdecvar('c',z);
Vdot = jacobian(V,x)*f;

sosc = polyconstr;
sosc(1) = s>=0;
sosc(2) = Vdot <= -1e-6*(x'*x) + s*(V+t);

% Solve with gsosopt
[info,dopt,sossol] = gsosopt(sosc,x,t,opts);
s = subs(s,dopt);
g = -info.tbnds(2); 
fprintf('gamma = %4.3f\n',g);

% Plot ROA estimate:= countour of V associated with gamma level set.
domain = [-3 3 -3 3];
title(['gamma = ' num2str(g)])
[C,ph2]=pcontour(V,g,domain,'r:');
set(ph2,'LineWidth',2)
hold off; grid on;

