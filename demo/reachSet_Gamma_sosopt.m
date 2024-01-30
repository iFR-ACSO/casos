close all
clear
clc


%% define variables
x = mpvar('x',2,1);
w = mpvar('w',1,1);
u = mpvar('u',1,1);

%% dynamics
f = [x(2);
     (1-x(1)^2)*x(2)-x(1)];

gd = [0;1];

xdot = f+gd*u;

%% Terminal Penalty & Constraints
% terminal LQR
[A,B] = plinearize(xdot, x, u);
[K,P] = lqr(A, B, eye(2), 2.5);

Pinit = lyap(A-B*K,eye(2));

%% constraints
umin = -1;
umax = 1;

% P = [6.4314 0.458; 0.458 5.8227];

l = x'*P*x-1;              % l(x) \leq 0
g = 3*x(2)^2 + x(1)^2 - 1;  % g(x) \leq 0

%% Problem Definitions
T     = 1;
t1    = 0;
t2    = T;

deg_V = 2;

maxIter = 100;


%% define SOS problem
pvar t w gamma

hT = (t-t1)*(t2-t);
hU = [(u-umin);(umax-u)];

s1  = sosdecvar('c1',monomials([t;x],0:2));
s2  = sosdecvar('c21',monomials([t;x],2));

s3 = sosdecvar('c3',monomials([t;x],0:2));
s4 = sosdecvar('c4',monomials([t;x],2));

s51  = sosdecvar('c51',monomials([t;x],0:2));
s52  = sosdecvar('c52',monomials([t;x],0:2));
s61  = sosdecvar('c61',monomials([t;x],2));
s62  = sosdecvar('c62',monomials([t;x],2));

s9   = sosdecvar('c9',monomials([t;x],2));
s10  = sosdecvar('c10',monomials([t;x],2));
s11  = sosdecvar('c11',monomials([t;x],2));

s7   = sosdecvar('c7' ,monomials(x,0:2));
s8   = sosdecvar('c8' ,monomials(x,0:2));

V = polydecvar('V',monomials([t;x],0:4));
k = polydecvar('k',monomials([t;x],0:4));

Vval = x'*P*x

sopt = sosoptions;
sopt.solver = 'mosek';
sopt.simplify = 'off';
% sopt.solveropts.param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean ups
% sopt.solveropts.param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

gsopt = gsosoptions;
gsopt.solver = 'mosek';
gsopt.maxobj = 0;
gsopt.minobj = -10;
gsopt.simplify = 'off';
gsopt.display ='on';
% gsopt.solveropts.param.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT'; % Use just the interior point algorithm to clean ups
% gsopt.solveropts.param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER'; % Don't use basis identification (it's slow)

% profile on -historysize 50000000000
gamma_old = [];

R = 0.5;
q = 1;
profile on -historysize 500000000000
tic
 for iter = 1:10
    % solve for multiplier and control law
    sosconk = [s1 >= 0;
               s2 >= 0; 
               s3 >= 0; 
               s4 >= 0; 
               s51 >= 0; 
               s52 >= 0; 
               s61 >= 0; 
               s62 >= 0;
               s8 >=  0;
               s8*(subs(Vval,t,t2)+gamma) - l >= 0;...
               subs(hU(1),u,k) + s51*(Vval+gamma) - s61*hT  >= 0;...
               subs(hU(2),u,k) + s52*(Vval+gamma) - s62*hT  >= 0;...
               s1*(Vval+gamma) - s2*hT - jacobian(Vval,t) - jacobian(Vval,x)*f - jacobian(Vval,x)*gd*k >= 0;...
               s3*(Vval+gamma) - s4*hT - g >= 0
               ];

[info,dopt] = gsosopt(sosconk,[x;t;w],gamma,gsopt);

if ~info.feas
    disp(['h step infeasible in iteration: ' num2str(iter)]); 
%     break
else
    disp(['h step feasible in iteration:' num2str(iter)])
    s1val  = subs(s1,dopt);
    s2val  = subs(s2,dopt);
    s3val  = subs(s3,dopt);
    s4val  = subs(s4,dopt);
    s51val = subs(s51,dopt);
    s52val = subs(s52,dopt);
    s61val = subs(s61,dopt);
    s62val = subs(s62,dopt);
    s8val  = subs(s8,dopt);
    kval   = subs(k,dopt);
    gammaVal = -double(subs(gamma,dopt));


end

 % solve for storage function
sosconV = [
           s7 >= 0;...
           s1val*(V-gammaVal) - s2*hT - jacobian(V-gammaVal,t) - jacobian(V-gammaVal,x)*f - jacobian(V-gammaVal,x)*gd*kval >= 0;...
           s3val*(V-gammaVal) - s4*hT - g >= 0;...
           subs(hU(1),u,kval) + s51val*(V-gammaVal) - s61*hT  >= 0;...
           subs(hU(2),u,kval) + s52val*(V-gammaVal) - s62*hT  >= 0;...
           s8val*(subs(V,t,t2)-gammaVal) - l >= 0
           s7*(subs(Vval,t,t1)-gammaVal)-(subs(V,t,t1)-gammaVal) >= 0
           ];


[info,dopt] = sosopt(sosconV,[x;t;w],sopt);
   
if ~info.feas
    disp(['V step infeasible in iteration:' num2str(iter)]); 
    break
else
    disp(['V step feasible in iteration:' num2str(iter)])
    Vval = subs(V,dopt);

end


% if isempty(gamma_old)
%     gamma_old = gammaVal;
% else
%     if abs(gamma_old-gammaVal) <= 1e-3
%         disp(['Gamma converged in iteration:' num2str(iter)])
%         break
%     else
%         gamma_old = gammaVal;
%     end

% end

 end
 toc
profile viewer
Vval = subs(Vval,t,t+t1)-gammaVal

Vval = cleanpoly(Vval,1e-4);

%% plotting
domain = [-2 2 -2 2];
figure()
pcontour(subs(Vval,t,0),0,domain,'b--')
hold on
pcontour(subs(Vval,t,T),0,domain,'b')
hold on
pcontour(l,0,domain,'r')
hold on
pcontour(g,0,domain,'k')
% domain = [x(1), -2, 2; x(2), -2, 2];
% pvolume(subs(Vval,t,0),0,domain)
% pvolume(subs(Vval,t,T),0,domain)
% legend('BRS','Inner Terminal','Terminal Set','State Constraint')

% for dt = 0.1:0.1:T-0.1
% pcontour(subs(Vval,t,dt),0,domain,'g--')   
% hold on
% end
