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
Dmax = diag([20 20*d2r 50*d2r 20*d2r]);

% target function
l = x'*Dmax'*blkdiag(1/(4)^2, 1/(pi/30)^2, 1/(pi/15)^2, 1/(pi/30)^2)*Dmax*x - 1;


Apre = nabla(x_dot,x);
Bpre = nabla(x_dot,u);

% equilibriums
xbar = zeros(4,1);
ubar = 0;

% substitute in the value of equilibrium
A = double(subs(subs(Apre,x,xbar),u,ubar));
B = double(subs(subs(Bpre,x,xbar),u,ubar));


Q = diag([1 1 1 1]);
R = diag(1);
K = lqr(A,B,Q,R);

P = lyap((A-B*K)',eye(4));

Vval = x'*P*x;
% 
% load('step3.mat','Vval_list')
% Vval = Vval_list(end);

% Vval = Vval.coefficient'*monomials([x;t],0:4)


% Trim point for elevator channel is 0.0489 rad.
% Saturation limit for elevator channel is -10 deg to 10 deg
um =  -(10*d2r + 0.0489); 
uM =    10*d2r - 0.0489;

% affine control constraints
u_min = u - um; 
u_max = uM - u;


%% define SOS problem polynomials

% reachability storage function
V = casos.PS.sym('v',monomials([x;t],0:4));
k = casos.PS.sym('k',monomials([x;t],0:4),[length(u) 1]);

% SOS multiplier dissipation inequality
s1 = casos.PS.sym('s1',monomials([x;t],0:4));
s2 = casos.PS.sym('s2',monomials([x;t],0:4));

% SOS multiplier for control constraints
s51 = casos.PS.sym('s51',monomials([x;t],0:4),[length(u) 1]);
s52 = casos.PS.sym('s52',monomials([x;t],0:4),[length(u) 1]);

s61 = casos.PS.sym('s61',monomials([x;t],0:4),[length(u) 1]);
s62 = casos.PS.sym('s62',monomials([x;t],0:4),[length(u) 1]);

% SOS multiplier terminal-set inclusion
s7 = casos.PS.sym('s7',monomials(x,0:4));  

% % SOS multiplier grow constraint
% s8 = casos.PS.sym('s8',monomials(x,0:2));   


g = casos.PS.sym('g');

% Time horizon
T  = 3;

t0 = 0;
t1 = T;

% time polynomial
hT = (t)*(T-t);


%% setup solver
disp('==================================================================')
disp('Setting up solver ... ')
tic

V_sym = casos.PS.sym('V',basis(V));

sos1 = struct('x',[k;s1;s51;s52;s2;s61;s62;s7],...
              'p',[V_sym;g]);     

sos1.('g') = [s1; s2;s51;s52;s61;s62;
              s7-0.0001;
              s1*(V_sym-g) - s2*hT - nabla((V_sym),t) - nabla((V_sym),x)*(f+gx*k);
              subs(u_min,u,k)     + s51.*(V_sym-g) - s61.*hT;
              subs(u_max,u,k)     + s52.*(V_sym-g) - s62.*hT
              subs(V_sym,t,T) - g  - s7*l
              ];

% states + constraint are SOS cones
% opts = struct('sossol','mosek');
% opts.Kx = struct('s', size(sos1.x( (length(u)+1):end) ,1), ...
%                  'l',length(u));                

opts.Kx = struct('s', 0, ...
                 'l',size(sos1.x,1));                


opts.Kc = struct('s', length(sos1.g));          % we only have sos con.

% setupsolve
opts.error_on_fail = 0; 

solve_Sstep = casos.sossol('S1','mosek',sos1,opts);

tend1 = toc;

disp(['setting up solver 1 took ' num2str(tend1) ' s'])

% % copy SOS constraints; currently SOS poly cannot be used for
% % parameterization; and helps code readability
% s1_sym = casos.PS.sym('s1',basis(s1));
% s2_sym = casos.PS.sym('s2',basis(s2));
% 
% s51_sym = casos.PS.sym('s51',basis(s51(1)),[length(u) 1]);
% s52_sym = casos.PS.sym('s52',basis(s52(1)),[length(u) 1]);
% 
% s61_sym = casos.PS.sym('s61',basis(s61(1)),[length(u) 1]);
% s62_sym = casos.PS.sym('s62',basis(s62(1)),[length(u) 1]);
% 
% s7_sym  = casos.PS.sym('s7',basis(s7));
% s8_sym  = casos.PS.sym('s8',basis(s8));
% 
% k_sym = casos.PS.sym('V',basis(k));
% tic
% sos2 = struct('x',[V;s2;s61;s62;s7;s8], ...
%               'p',[g;V_sym;k_sym;s1_sym;s51_sym;s52_sym]);     
% 
% sos2.('g') = [
%               s1_sym*(V-g) - s2*hT - nabla(V,t) - nabla(V,x)*(f+g*k_sym);
%               subs(u_min,u,k_sym)  + s51_sym.*(V-g) - s61.*hT;
%               subs(u_max,u,k_sym)  + s52_sym.*(V-g) - s62.*hT;
%               (s7-0.0001)*l - (subs(V,t,t1)-g); 
%               s8*(subs(V_sym,t,t0)-g)  + g - subs(V,t,t0);
%               ];

% 
% % states + constraint are SOS cones
% opts    = struct;
% opts.Kx = struct('s', length(sos2.x(2:end)), ... % first entry storage fun.
%                  'l',1); 
% opts.Kc = struct('s', length(sos2.g));           % we only have sos con.
% 
% % setup solver
% solver_Vstep = casos.sossol('S2','mosek',sos2,opts);

% tend2 = toc;



%% V-s-iteration
% disp(['setting up solver 2 took ' num2str(tend2) ' s'])

disp('==================================================================')

tic

itermax = 10;

disp('Start V-S-iteration ...')
for iter = 1:itermax

   gamma_ub = 1;
    gamma_lb = 0;
    num_exp = 0;
    go = 1;
    while (num_exp <= 12 || go)
        num_exp = num_exp + 1;
        if num_exp >= 20
            break
        end
        
        gamma_try = (gamma_ub + gamma_lb)/2

        % find mulitplier and control law
        sol_Sstep = solve_Sstep('p',[Vval;gamma_try]);
    
        switch (solve_Sstep.stats.UNIFIED_RETURN_STATUS)
            case 'SOLVER_RET_SUCCESS' 
                 gamma_lb = gamma_try;
                go = 0;
                gamma_use = gamma_try;
                disp(['S-step feasible in ' num2str(iter) '/' num2str(itermax)])
                
            case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
    
                disp(['g-step infeasible in ' num2str(iter) '/' num2str(itermax)])
                   gamma_ub = gamma_try;
        end
    end

   % find reachability storage function
   sol_Vstep = solver_Vstep('p',[-sol_Sstep.f;Vval;sol_Sstep.x(1:7)]);


   switch (solver_Vstep.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            Vval = sol_Vstep.x(1);

            disp(['V-step feasible in ' num2str(iter) '/' num2str(itermax) ])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter) '/' num2str(itermax)])
            break
                    
        otherwise, error('Failed.')
   end

end

tend4  = toc;
disp(['Finished V-s-iteration after ' num2str(tend4) ])




