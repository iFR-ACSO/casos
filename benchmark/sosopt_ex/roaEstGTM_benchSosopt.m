function [gval,bval,solveTime,solverTime,buildTime]= roaEstGTM_benchSosopt()

%System dynamics
pvar x1 x2 x3 x4

x= [x1; x2; x3;x4];


%% Polynomial Dynamics
f1 = @(V,alpha,q,theta,eta,F) +1.233e-8*V.^4.*q.^2         + 4.853e-9.*alpha.^3.*F.^3    + 3.705e-5*V.^3.*alpha.*q      ...
                             - 2.184e-6*V.^3.*q.^2         + 2.203e-2*V.^2.*alpha.^3     - 2.836e-6*alpha.^3.*F.^2      ...
                             + 3.885e-7*alpha.^2.*F.^3     - 1.069e-6*V.^3.*q            - 4.517e-2*V.^2.*alpha.^2      ...
                             - 2.140e-3*V.^2.*alpha.*eta   - 3.282e-3*V.^2.*alpha.*q     - 8.901e-4*V.^2.*eta.^2        ...
                             + 9.677e-5*V.^2.*q.^2         - 2.037e-4*alpha.^3.*F        - 2.270e-4*alpha.^2.*F.^2      ...
                             - 2.912e-8*alpha.*F.^3        + 1.591e-3*V.^2.*alpha        - 4.077e-4*V.^2.*eta           ...
                             + 9.475e-5*V.^2.*q            - 1.637   *alpha.^3           - 1.631e-2*alpha.^2.*F         ...
                             + 4.903   *alpha.^2.*theta    - 4.903   *alpha.*theta.^2    + 1.702e-5*alpha.*F.^2         ...
                             - 7.771e-7*F.^3               + 1.634   *theta.^3           - 4.319e-4*V.^2                ...
                             - 2.142e-1*alpha.^2           + 1.222e-3*alpha.*F           + 4.541e-4*F.^2                ...
                             + 9.823   *alpha              + 3.261e-2*F                  - 9.807.*theta      + 4.284e-1;
% change of pitch angle
f4 = @(V,alpha,q,theta,eta,F)  q;
                      
% change of angle of attack
f2 = @(V,alpha,q,theta,eta,F) -3.709e-11*V.^5.*   q.^2     + 6.869e-11*V.*alpha.^3.*F.^3 + 7.957e-10*V.^4.*alpha.*q     ...
                             + 9.860e-9 *V.^4.*   q.^2     + 1.694e-5 *V.^3.*alpha.^3    - 4.015e-8 *V.*alpha.^3.*F.^2  ...
                             - 7.722e-12*V.*alpha.^2.*F.^3 - 6.086e-9 *alpha.^3.*F.^3    - 2.013e-8 *V.^4.*q            ...
                             - 5.180e-5 *V.^3.*alpha.^2    - 2.720e-6 *V.^3.*alpha.*eta  - 1.410e-7 *V.^3.*alpha.*q     ...
                             + 7.352e-7 *V.^3.*eta.^2      - 8.736e-7 *V.^3.*q.^2        - 1.501e-3 *V.^2.*alpha.^3     ...
                             - 2.883e-6 *V.*alpha.^3.*F    + 4.513e-9 *V.*alpha.^2.*F.^2 - 4.121e-10*V.*alpha.*F.^3     ...
                             + 3.557e-6 *alpha.^3.*F.^2    + 6.841e-10*alpha.^2.*F.^3    + 4.151e-5 *V.^3.*alpha        ...
                             + 3.648e-6 *V.^3.*eta         + 3.566e-6 *V.^3.*q           + 6.246e-6 *V.^2.*alpha.*q     ...
                             + 4.589e-3 *V.^2.*alpha.^2                                  - 6.514e-5 *V.^2.*eta.^2       ...
                             + 2.580e-5 *V.^2.*q.^2        - 3.787e-5 *V.*alpha.^3       + 3.241e-7 *V.*alpha.^2.*F     ...
                             + 2.409e-7 *V.*alpha.*F.^2    + 1.544e-11*V.*F.^3           + 2.554e-4 *alpha.^3.*F        ...
                             - 3.998e-7 *alpha.^2.*F.^2    + 3.651e-8 *alpha.*F.^3       + 4.716e-7 *V.^3               ...
                             - 3.677e-3 *V.^2*alpha        - 3.231e-4 *V.^2.*eta         - 1.579e-4 *V.^2.*q            ...
                             + 2.605e-3 *V.*alpha.^2       + 1.730e-5 *V.*alpha.*F       - 5.201e-3 *V.*alpha.*theta    ...
                             - 9.026e-9 *V.*F.^2           + 2.601e-3 *V.*theta.^2       + 3.355e-3 *alpha.^3           ...
                             - 2.872e-5 *alpha.^2.*F       - 2.134e-5 *alpha.*F.^2       - 1.368e-9 *F.^3               ...
                             - 4.178e-5 *V.^2              + 2.272e-4 *V.*alpha          - 6.483e-7 *V.*F               ...
                             - 2.308e-1 *alpha.^2          - 1.532e-3 *alpha.*F          + 4.608e-1 *alpha.*theta       ...
                             - 2.304e-1 *theta.^2          + 7.997e-7 *F.^2              - 5.210e-3 *V                  ...
                             - 2.013e-2 *alpha             + 5.744e-5 *F          + q    + 4.616e-1;
                         
% change of pitch rate                         
f3 = @(V,alpha,q,theta,eta,F) -6.573e-9*V.^5.*    q.^3     + 1.747e-6*V.^4.*q.^3         - 1.548e-4*V.^3.*q.^3          ...
                             - 3.569e-3*V.^2.*alpha.^3     + 4.571e-3*V.^2.*q.^3         + 4.953e-5*V.^3.*q             ...
                             + 9.596e-3*V.^2.*alpha.^2     + 2.049e-2*V.^2.*alpha.*eta   - 2.431e-2*V.^2.*alpha         ...
                             - 3.063e-2*V.^2.*eta          - 4.388e-3*V.^2.*q            - 2.594e-7*F.^3                ...
                             + 2.461e-3*V.^2               + 1.516e-4*F.^2               + 1.089e-2*F        + 1.430e-1;                         

Kq = 0.0698;

%% Trim condition
v0      = 45;       % m/s
alpha0  = 0.04924; % rad
q0      = 0;       % rad/s
theta0  = 0.04924; % rad

eta0    = 0.04892; % rad
delta0  = 14.33;   % percent
            

x0 = [v0; alpha0; q0; theta0];
u0 = [eta0; delta0];

% get trim point 
f = @(x,u) [
                    f1(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
                    f2(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
                    f3(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
                    f4(x(1,:),x(2,:),x(3,:),x(4,:),u(1,:),u(2,:))
];
[x0, u0] = findtrim(f,x0, u0);


% set up dynamic system
f = f(x+x0,[Kq*x(3); 0]+u0);

d = ([
          convvel(20, 'm/s', 'm/s')  %range.tas.lebesgue.get('ft/s')
          convang(20, 'deg', 'rad')  %range.gamma.lebesgue.get('rad')
          convang(50, 'deg', 'rad')  %range.qhat.lebesgue.get('rad')
          convang(20, 'deg', 'rad')  %range.alpha.lebesgue.get('rad')
]);


D = diag(d)^-1;
f = subs(D*f, x, D^-1*x);

f = cleanpoly(f, 1e-6, 0:5);

A = jacobian(f,x);

A0 = double(subs(A,x,zeros(4,1))); 

% solve Lyapunov equation
P = lyap(A0',10*eye(4));

Vval = x'*P*x;


p = x'*x*1e2;

% epsilon polynomial
l = x'*x*1e-6;

bval = [];
gval = [];

endTimeParse1 = [];
endTimeParse2 = [];

solverTime1 = [];
solverTime2 = [];
solverTime3 = [];


% Form bilinear constraints.  The variable t:=-gamma is used to convert
% the maximization of gamma into a minimization of t.
pvar t b;
V  = polydecvar('v',monomials(x, 2:4 ));
s1 = sosdecvar('s1',monomials(x, 1:2 ));
s2 = sosdecvar('s2',monomials(x, 0:2 ));

bval_old = [];

solveTime_start = tic;
for iter = 1:100

opts = gsosoptions;
opts.minobj = -1000; 
opts.maxobj = 0;
% opts.display ='on';
% opts.simplify = 'off';
opts.solver = 'mosek';
% opts.simplify = 'off';
% opts.scaling = 'on';

    % gamma-step    
    startTimeParse1 = tic;
    sosc = polyconstr;
    sosc(1) = s1 >=0;
    sosc(2) = s1*(Vval + t) - jacobian(Vval,x)*f - l >= 0;
    
    % Solve with gsosopt
    [info,dopt] = gsosopt(sosc,x,t,opts);

    s1val = subs(s1,dopt);
    gval = -info.tbnds(2); 
    
    solverTime1(iter) = sum(arrayfun(@(x) x.solverinfo.solvertime, info.sdpsol));
    endTimeParse1     = [endTimeParse1 toc(startTimeParse1)-solverTime1(iter)];


     % beta-step    
    startTimeParse2 = tic;
    sosc = polyconstr;
    sosc(1) = s2 >=0;
    sosc(2) = s2*(p + b) + gval - Vval >= 0;
    
    % Solve with gsosopt
    [info,dopt] = gsosopt(sosc,x,b,opts);
    % 
    s2val = subs(s2,dopt);
    bval = -info.tbnds(2); 
    
    solverTime2(iter) = sum(arrayfun(@(x) x.solverinfo.solvertime, info.sdpsol));
    endTimeParse2    = [endTimeParse2 toc(startTimeParse2)-solverTime2(iter)];


     % beta-step    
    startTimeParse3 = tic;
    sosc = polyconstr;
    sosc(1) = V - l >=0;
    sosc(2) = s1val*(V - gval) - jacobian(V,x)*f - l >= 0;
    sosc(3) = s2val*(p - bval) + gval - V >= 0;
    
    % Solve with gsosopt
    opts = sosoptions;
    opts.solver = 'mosek';
    % opts.simplify = 'off';

    [info,dopt,~] = sosopt(sosc,x,opts);
    solverTime3(iter) = info.sdpsol.solverinfo.solvertime;
    endTimeParse3(iter) =  toc(startTimeParse3)-solverTime3(iter);
    Vval = subs(V,dopt);


    fprintf('Iteration %d: b = %g, g = %g.\n',iter,full(bval),full(gval));
	
				
  
   if ~isempty(bval_old)
        if abs(full(bval-bval_old)) <= 1e-3
            break
        else
            bval_old = bval;
        end
    else
        bval_old = bval;
    end
end

solveTime  = toc(solveTime_start);
buildTime  = sum(endTimeParse1) + sum(endTimeParse2) + sum(endTimeParse3);
solverTime = sum(solverTime1) + sum(solverTime2) + sum(solverTime3);


