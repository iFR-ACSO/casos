%--------------------------------------------------------------------------
% 
% Implementation of custom V-s-iteration for the GTM 4D ROA problem in 
% CaSoS. The quasi-convex solver is used to perform the bisections.
%
%--------------------------------------------------------------------------

function [gval, solverTime_total,buildTime] = reachEstGTM_benchCasos()

% system states 
x = casos.Indeterminates('x',4);

x1 = x(1); % V 
x2 = x(2); % alpha
x3 = x(3); % q
x4 = x(4); % theta
t = casos.Indeterminates('t');

% Polynomial Dynamics
d2r = pi/180;
% scale matrix
Dmax = diag([20 20*d2r 50*d2r 20*d2r]);

% scaled 4-state system
f = [ % f1
      - 0.01955829859645207*x1^2*x2 + 0.0006116050110168922*x1^2*x3 - 0.4597074905072323*x1*x2^2 ...
      - 0.02143363124979007*x1*x2*x3 + 0.0913633506074555*x2^3 + 0.0104276150041391*x2^2*x4 ...
      - 0.0104276150041391*x2*x4^2 + 0.003475871668046367*x4^3 - 0.00975287614756193*x1^2 ...
      - 0.08801234368403432*x1*x2 - 0.001251382140034897*x1*x3 - 0.5201826057909725*x2^2 ...
      - 0.04982793763401992*x2*x3 + 0.00180695498489715*x3^2 - 0.04388794266402869*x1 ...
      + 0.07194181873717953*x2 - 0.00290915666614335*x3 - 0.1711592037553279*x4;...
      % f2
      0.05879789171572614*x1^3 + 0.6755306062430397*x1^2*x2 + 0.07878176650294093*x1^2*x3 ...
      + 0.6715746030304725*x1*x2^2 - 0.03631259653379534*x1*x2*x4 - 0.006619604539373272*x1*x3^2 ...
      + 0.01815629826689767*x1*x4^2 - 0.1817374041405456*x2^3  + 0.1364168103441202*x1^2 ...
      - 1.371702883566865*x1*x2 + 0.005525130555690897*x1*x3 + 1.479555253341581*x2^2 ...
      + 0.002672981895904953*x2*x3 + 0.07915793925313672*x2*x4 + 0.0135837161533212*x3^2 ...
      - 0.03957896962656836*x4^2  - 0.5767816606949114*x1 - 3.236059303865032*x2 + 2.30669948822443*x3;...
      % f3
      - 3.582346855286648*x1^2*x2 + 0.9194217971573437*x1^2*x3 + 2.279359758657476*x1*x2^2 ...
      - 0.3522980647489417*x2^3 - 0.07456233935706458*x1^2 - 16.12056084878992*x1*x2 ...
      - 1.881194554322757*x1*x3 + 2.56427972848966*x2^2  - 0.33553052710679*x1 ...
      - 18.13563095488866*x2 - 4.373316114187153*x3;...
      % f4
      2.5*x3];
g = [% g1
     - 0.01200284890938563*x1^2- 0.0672488718061946*x1*x2- 0.05401282009223534*x1 ...
     - 0.07565498078196893*x2- 0.06076442260376476; ...
     % g2
     0.190652302942119*x1^2 + 0.1033082169845372*x1*x2 - 0.3900865469718997*x1 ...
     + 0.240166275746567*x2 - 0.9068555816726751; ...
     % g3
     -13.57875281932966*x1^2 + 14.75731515054161*x1*x2 - 61.10438768698349*x1 ...
     + 16.60197954435931*x2 - 68.74243614785642; ...
     % g4
     0
    ];

% target function
rT = x'*Dmax'*blkdiag(1/(4)^2, 1/(pi/30)^2, 1/(pi/15)^2, 1/(pi/30)^2)*Dmax*x - 1;


P = [4.54064415441163	0.139789072082587	-0.0107315729582360	-0.406387155887787
     0.139789072082587	0.165595957230149	-0.00115647767230559	-0.0841739454608939
    -0.0107315729582360	-0.00115647767230559	0.00732884319609208	0.00922486739143782
   -0.406387155887787	-0.0841739454608939	0.00922486739143782	0.514580693591350];

Vval = x'*P*x;

% Trim point for elevator channel is 0.0489 rad.
% Saturation limit for elevator channel is -10 deg to 10 deg
uM = 10*d2r - 0.0489;
um = -(10*d2r + 0.0489);

% start to measure build time of all parameterized solver
buildTime_start = tic;

% Lyapunov function candidate
V    = casos.PS.sym('v',monomials([x;t],0:4));
Vold = casos.PS.sym('vo',sparsity(V));
K    = casos.PS.sym('k',monomials([x;t],0:4));

% SOS multiplier
s1 = casos.PS.sym('s1',monomials(x,0:2),'gram');
s2 = casos.PS.sym('s2',monomials([x;t],0:2),'gram');
s3 = casos.PS.sym('s3',monomials([x;t],0:2),'gram');
s4 = casos.PS.sym('s4',monomials(x,0:2),'gram');
s5 = casos.PS.sym('s5',monomials([x;t],0:2),'gram');
s6 = casos.PS.sym('s6',monomials([x;t],0:2),'gram');
s7 = casos.PS.sym('s7',monomials([x;t],0:2),'gram');
s8 = casos.PS.sym('s8',monomials([x;t],0:2),'gram');

% level of stability
b = casos.PS.sym('b');

T = 3;
h = casos.PS(1*t*(T-t*1));

% options
opts               = struct('sossol','mosek');
opts.error_on_fail = 0;
opts.conf_interval = [-1 0];

%% setup solver
sos1 = struct('x',[K;s2;s3;s4;s5;s6;s7;s8], ... % dec.var
              'f',-b, ... % cost function for bisection
              'p',V);     % parameter

% constraint
sos1.('g') = [s4-0.0001;...
             -(nabla(V, t) + nabla(V, x)*(f + g*K)) - s2*h + s3*(V - b); ...
             -s4*rT + subs(V,t,T) - b;...
             uM - K + s5*(V - b) - s6*h; ...
             K - um + s7*(V - b) - s8*h];


% states + constraint are SOS cones
opts.Kx = struct('lin', 1,'sos',7);
opts.Kc = struct('sos', 5);

% build first solver
S1 = casos.qcsossol('S1','bisection',sos1,opts);


% solver 2: V-step
sos2 = struct('x',[V;s1;s2;s4;s6;s8], ...        % dec.var
              'p',[b;K;s3;s5;s7;Vold]); % parameter

% constraints
sos2.('g') = [  s4-0.0001;...
                -(subs(V,t,0)-b) + s1*(subs(Vold,t,0)-b);
                -(nabla(V, t) + nabla(V, x)*f + K*nabla(V, x)*g) - s2*h + s3*(V - b);...
                -s4*rT + subs(V,t,T) - b;...
                uM - K + s5*(V - b) - s6*h; ...
                K - um + s7*(V - b) - s8*h];

opts    = struct;
opts.Kx = struct('sos', 5, 'lin', 1); 
opts.Kc = struct('sos', 6);

% build second solver
S2 = casos.sossol('S','mosek',sos2,opts);

tmpbuildTime = toc(buildTime_start);

% initialize arrays
solvetime_all1 = zeros(100,1);
solvetime_all2 = zeros(100,1);
solvetime_all3 = zeros(100,1);


buildTime1 = zeros(100,1);
buildTime2 = zeros(100,1);

%% V-s-iteration
for iter = 1:10

    % gamma step
     startSolve1 =tic;
    sol1 = S1('p',Vval);    
    % get all solver times from all subiterations of the biscetion
    solvetime_all1(iter) = sum(cellfun(@(x) x.solvetime_matlab, S1.stats.iter));
    buildTime1(iter) = toc(startSolve1) -solvetime_all1(iter);

    % extract solution
    if strcmp(S1.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') 
        gval = -sol1.f;
        Kval =  sol1.x(1);
        s3val = sol1.x(3);
        s5val = sol1.x(5);
        s7val = sol1.x(7);
    else
        disp('Problem is infeasible in gamma-step! ')
        break
    end

    % beta step
    startSolve2 =tic;
    sol2 = S2('p',[gval;Kval;s3val;s5val;s7val;Vval]); 
    % get solver time
    solvetime_all2(iter) = S2.stats.solvetime_matlab;
    buildTime2(iter) = toc(startSolve2) -solvetime_all2(iter);
        
    
     if strcmp(S2.stats.UNIFIED_RETURN_STATUS,'SOLVER_RET_SUCCESS') 
           % extract solution
            Vval = sol2.x(1);
     else
           disp('Problem is infeasible in V-step! ')
        break
     end

    
    % show progress 
    fprintf('Iteration %d: , g = %g.\n',iter,full(gval));
    
    % % check convergence
    % if ~isempty(bval_old)
    %     if abs(full(bval-bval_old)) <= 1e-3
    %         break
    %     else
    %         bval_old = bval;
    %     end
    % else
    %     bval_old = bval;
    % end

end % end for-loop

% total solver time over all iterations
buildTime  = tmpbuildTime + sum(buildTime1) + + sum(buildTime2);
solverTime_total = sum(solvetime_all1) + sum(solvetime_all2) + sum(solvetime_all3);

end % end of function