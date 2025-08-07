%% ------------------------------------------------------------------------
%
%
%   Short Descirption:  Calculate an inner-estimate of the
%                       backward reachable set for the longitudinal motion 
%                       of the Nasa Generic Transport Model. To increase
%                       the size of the sublevel set we try to minimize the
%                       squared distance to a defined set. Scaled model,
%                       constraints and time horizon from [1] and [2].
%
%   References: 
%          [1] Yin, H., Arcak, M. and Seiler, P. ,
%              Backward Reachability for Polynomial Systems on a Finite Horizon, 
%              IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 66, NO. 12, DECEMBER 2021,
%              doi: 10.1109/TAC.2021.3056611
%          [2] https://github.com/heyinUCB/Backward-Reachability-Analysis-and-Control-Synthesis 
%
%--------------------------------------------------------------------------


clear
clc

% system states 
x = casos.Indeterminates('x',4);
t = casos.Indeterminates('t');

x1 = x(1); % V 
x2 = x(2); % alpha
x3 = x(3); % q
x4 = x(4); % theta

% Polynomial Dynamics
d2r = pi/180;
% scale matrix see reference [2]
Dmax = diag([20 20*d2r 50*d2r 20*d2r]);

% scaled 4-state system see reference [2]
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

% terminal set see reference [2]
rT = x'*Dmax'*blkdiag(1/(4)^2, 1/(pi/30)^2, 1/(pi/15)^2, 1/(pi/30)^2)*Dmax*x - 1;


P = [4.54064415441163	0.139789072082587	-0.0107315729582360	-0.406387155887787
     0.139789072082587	0.165595957230149	-0.00115647767230559	-0.0841739454608939
    -0.0107315729582360	-0.00115647767230559	0.00732884319609208	0.00922486739143782
   -0.406387155887787	-0.0841739454608939	0.00922486739143782	0.514580693591350];

g0 = x'*P*x-1;

% Trim  for elevator and control limits see reference [2]
uM =   10*d2r - 0.0489;
um = -(10*d2r + 0.0489);


% Lyapunov function candidate
V    = casos.PS.sym('v', monomials([x;t],0:4));  %
K    = casos.PS.sym('k',monomials([x;t],1));

% SOS multiplier
s2 = casos.PS.sym('s2',monomials([x;t],0:2));
s3 = casos.PS.sym('s3',monomials([x;t],0:4));
s4 = casos.PS.sym('s4',monomials([x;t],0:2));
s5 = casos.PS.sym('s5',monomials([x;t],0:2));
s6 = casos.PS.sym('s6',monomials([x;t],0:2));
s7 = casos.PS.sym('s7',monomials([x;t],0:2));
s8 = casos.PS.sym('s8',monomials([x;t],0:2));
s10 = casos.PS.sym('s10',monomials([x;t],0:2));
s9 = casos.PS.sym('s9',monomials(x,0:2));
b = casos.PS.sym('b');

% desired time horizon and time polynomial
T = 3;
h = casos.PS(1*t*(T-t*1));


%% setup solver
sos1 = struct('x',[K;s2;s3;s4;s5;s6;s7;s8;s9;s10], ...  % dec.var
              'p',[V;b]);                                  % parameter

% constraint
sos1.('g') = [
                s2;
                s3;
                s4;
                s5;
                s7;
                s6;
                s8;
                s9;
                s10;
                -(nabla(V, t) + nabla(V, x)*(f + g*K)) - s2*h + s3*(V-b); ...
                s4*(V-b) - s10*h - g0 ;...
                s9*(subs(V,t,T)-b) - subs(V,t,0) + b;...
                uM - K + s5*(V - b)   - s6*h; ...
                K - um + s7*(V -b )  - s8*h
             ];


% states + constraint are SOS cones
opts.Kx = struct('lin',length(sos1.x));
opts.Kc = struct('sos', length(sos1.g));

% build first solver
startBuildS1 = tic;
S1 = casos.sossol('S1','mosek',sos1,opts);
buildTime1 = toc(startBuildS1);


sos2 = struct('x',[V;b], ...                           % dec.var
              'f',dot(g0-(V-b),g0-(V-b)), ...      % cost function 
              'p',[K;s2;s3;s4;s5;s6;s7;s8;s9;s10]);    % parameter

% constraint
sos2.('g') = [ -(nabla(V, t) + nabla(V, x)*(f + g*K)) - s2*h + s3*(V-b); ...
                s4*(V-b) - s10*h - g0 ;...
                s9*(subs(V,t,T)-b) - subs(V,t,0) + b;...
                uM - K + s5*(V - b)   - s6*h; ...
                K - um + s7*(V - b )  - s8*h
             ];


% states + constraint are SOS cones
opts.Kx = struct('lin',length(sos2.x));
opts.Kc = struct('sos', length(sos2.g));
startBuildS2 = tic;
% build first solver
S2 = casos.sossol('S2','mosek',sos2,opts);
buildTime2 = toc(startBuildS2);

% initial storage function
Vval = x'*P*x;

% fix level set; for b >= 0.01 first step is infeasible
bval = 0.001;

fval_old = [];

startSolve = tic;
for iter = 1:100

% call solver 1 i.e. compute multiplier and control law
sol1 = S1('p',[Vval;bval]);  

% bval = 0.1;
sol2 = S2('p',[sol1.x]);

Vval = sol2.x(1);
bval = sol2.x(2);

% show progress 
fprintf('Iteration %d: f = %g, b = %g.\n',iter,full(sol2.f),full(bval));

% check convergence of cost function
 if ~isempty(fval_old)
    if abs(full(sol2.f-fval_old)) <= 1e-4
        
        break
    else
        fval_old = sol2.f;
    end
else
    fval_old = sol2.f;
 end


end
buildtime = buildTime1 + buildTime2
totalSolveTime = toc(startSolve)
timePerIteration = totalSolveTime/iter

n_con1 = S1.stats.conic.size_A.size(1)
n_dec1 = S1.stats.conic.size_A.size(2)

n_con2 = S2.stats.conic.size_A.size(1)
n_dec2 = S2.stats.conic.size_A.size(2)


%% plotting

Vval0 = subs(sol2.x(1),t,0)-bval;
VvalT = subs(sol2.x(1),t,T)-bval;

import casos.toolboxes.sosopt.*

figure(1)
clf
% reachability storage function at t = 0s
pcontour(subs(Vval0,x(3:4),zeros(2,1)),0,[-1 1 -1 1]*2,'b')
hold on 
pcontour(subs(VvalT,x(3:4),zeros(2,1)),0,[-1 1 -1 1]*2,'r')
pcontour(subs(g0,x(3:4),zeros(2,1)),0,[-1 1 -1 1]*2,'k--')
legend('V(0,x)','V(T,x)','r_T')



