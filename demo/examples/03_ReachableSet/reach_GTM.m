% ------------------------------------------------------------------------
%
%   Short Description:  Inner-approximation of the reachable set of the GTM
%                       model with respect to control constraint and
%                       terminal set.
%
%   Source: Yin et al. 
%
%   Date: 02/17/2024
%
% ------------------------------------------------------------------------


clear
clc
close all


%% define variables
x = casos.PS('x',4,1);
u = casos.PS('u',1,1);

% time
t  = casos.PS('t');

x1 = x(1); % V 
x2 = x(2); % alpha
x3 = x(3); % q
x4 = x(4); % theta



%% define your system here


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
gx = [% g1
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

% dynamics
x_dot =  f+gx*u;

d2r = pi/180;

% scale matrix
% Dmax = diag([20 20*d2r 50*d2r 20*d2r]);

% target function
% l = x'*Dmax'*blkdiag(1/(4)^2, 1/(pi/30)^2, 1/(pi/15)^2, 1/(pi/30)^2)*Dmax*x - 1;


Apre = nabla(x_dot,x);
Bpre = nabla(x_dot,u);

% equilibriums
xbar = zeros(4,1);
ubar = 0;

% substitute in the value of equilibrium; scaled dynamics was used
A = double(subs(subs(Apre,x,xbar),u,ubar));
B = double(subs(subs(Bpre,x,xbar),u,ubar));


Q = diag([1 1 1 1]);
R = diag(1);
K = lqr(A,B,Q,R);

P = lyap((A-B*K)',eye(4));

Vval0 = x'*P*x;



% Trim point for elevator channel is 0.0489 rad.
% Saturation limit for elevator channel is -10 deg to 10 deg
um =  -(10*d2r + 0.0489); 
uM =    10*d2r - 0.0489;

% affine control constraints
u_min = u - um; 
u_max = uM - u;

%% est. ROA to get a valid terminal set
xdot_c = subs(x_dot,u,-K*x);

[Vval,gval] = estROAGTM(xdot_c,Vval0,x);

l = Vval-gval;
Vval = Vval;
%% define SOS problem polynomials

% reachability storage function
V = casos.PS.sym('v',monomials([x;t],2));
k = casos.PS.sym('k',monomials([x;t],0:2),[length(u) 1]);

% SOS multiplier dissipation inequality
s1 = casos.PS.sym('s1',monomials([x;t],0:2));
s2 = casos.PS.sym('s2',monomials([x;t],2));

% SOS multiplier for control constraints
s51 = casos.PS.sym('s51',monomials([x;t],0:2),[length(u) 1]);
s52 = casos.PS.sym('s52',monomials([x;t],0:2),[length(u) 1]);

s61 = casos.PS.sym('s61',monomials([x;t],0:2),[length(u) 1]);
s62 = casos.PS.sym('s62',monomials([x;t],0:2),[length(u) 1]);

g = casos.PS.sym('g');

% Time horizon
T  = 1;

t0 = 0;
t1 = T;

% time polynomial
hT = (t)*(T-t);


%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic

V_sym = casos.PS.sym('V',basis(V));

sos1 = struct('x',[k;s1;s2;s51;s52;s61;s62],...
              'p',V_sym);     

sos1.('g') = [s1;s2;
              s1*(V_sym-gval) - s2*hT - nabla((V_sym),t) - nabla((V_sym),x)*(f+gx*k);
              subs(u_min,u,k)     + s51.*(V_sym-gval) - s61.*hT;
              subs(u_max,u,k)     + s52.*(V_sym-gval) - s62.*hT
              ];

opts.Kx = struct('s', 0, ...
                 'l',size(sos1.x,1));                


opts.Kc = struct('s', length(sos1.g));          % we only have sos con.

% setupsolve
opts.error_on_fail = 0; 

solve_Sstep = casos.sossol('S1','mosek',sos1,opts);

tend1 = toc;

disp(['setting up solver 1 took ' num2str(tend1) ' s'])


disp('==================================================================')

tic

itermax = 10;

disp('Start V-S-iteration ...')
for iter = 1:itermax


        % find mulitplier and control law
        sol_Sstep = solve_Sstep('p',Vval);
    
        switch (solve_Sstep.stats.UNIFIED_RETURN_STATUS)
            case 'SOLVER_RET_SUCCESS' 
             
                disp(['S-step feasible in ' num2str(iter) '/' num2str(itermax)])
                
            case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
    
                disp(['g-step infeasible in ' num2str(iter) '/' num2str(itermax)])
         
        end
   

   % % find reachability storage function
   % sol_Vstep = solver_Vstep('p',[-sol_Sstep.f;Vval;sol_Sstep.x(1:7)]);
   % 
   % 
   % switch (solver_Vstep.stats.UNIFIED_RETURN_STATUS)
   %      case 'SOLVER_RET_SUCCESS' 
   % 
   %          Vval = sol_Vstep.x(1);
   % 
   %          disp(['V-step feasible in ' num2str(iter) '/' num2str(itermax) ])
   % 
   %      case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
   % 
   %          disp(['V-step infeasible in ' num2str(iter) '/' num2str(itermax)])
   %          break
   % 
   %      otherwise, error('Failed.')
   % end

end

tend4  = toc;
disp(['Finished V-s-iteration after ' num2str(tend4) ])




