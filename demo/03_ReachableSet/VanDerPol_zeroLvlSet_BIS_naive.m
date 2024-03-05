% Inner-approximation reachable set Van-der-Pol Oscillator
close all
clear
clc

% system states
x = casos.PS('x',2,1);

t  = casos.PS('t');
t0 = casos.PS.sym('t0');
t1 = casos.PS.sym('t1');

% system dynamics
f = [x(2);
        (1-x(1)^2)*x(2)-x(1)];

gx = [0;1];


% terminal region
P = [6.4314    0.4580
    0.4580    5.8227];

l = x'*P*x-1;              % l(x) \leq 0

% state constraint
g = 3*x(2)^2 + x(1)^2 -1;  % g(x) \leq 0

% control constraint
umin = -1;
umax = 1;

% Time horizon
T  = 1;
dt = 0.1;

N = 100;

maxIter = N;

% time polynomial
hT = (t)*(dt-t);

% Storage function and control law
V = casos.PS.sym('v',monomials([x;t],0:6));
k = casos.PS.sym('k',monomials(x,0:4));

lsym  = casos.PS.sym('l',basis(subs(V,t,t1)));

% SOS multiplier

s1 = casos.PS.sym('s1',monomials([x;t],0:4),'gram'); % dissipation ineq.
s2 = casos.PS.sym('s2',monomials([x;t],4),'gram');  % dissipation ineq.

s3 = casos.PS.sym('s3',monomials([x;t],0:2),'gram');% state constraint
s4 = casos.PS.sym('s4',monomials([x;t],2),'gram');  % state constraint

s51 = casos.PS.sym('s51',monomials([x;t],2),'gram'); % control constraint
s52 = casos.PS.sym('s52',monomials([x;t],2),'gram'); % control constraint

s61 = casos.PS.sym('s61',monomials([x;t],2),'gram');   % control constraint
s62 = casos.PS.sym('s62',monomials([x;t],2),'gram');   % control constraint

s7 = casos.PS.sym('s7',monomials(x,0:1),'gram');   % grow constraint
s8 = casos.PS.sym('s8',monomials(x,0:1),'gram');   % grow constraint


%% setup solver

% solver 1: sos-multiplier and control law
tic
disp('-----------------------------------------')
disp('Building solver ... ')

sos1 = struct('x',[k;s1;s2;s3;s4;s51;s52;s61;s62],'p',[V;t0;t1]);     % parameter

sos1.('g') = [s1*V - s2*hT - nabla(V,t) - nabla(V,x)*(f + gx*k);
              k - umin  + s51*V - s61*hT;
              umax - k   + s52*V - s62*hT;
              s3*V - s4*hT - g];

% states + constraint are SOS cones
opts    = struct;
opts.Kx = struct('s', 8,'l',1);
opts.Kc = struct('s', 4);

klb = casos.PS(basis(k),-inf);
kub = casos.PS(basis(k),+inf);

% solver for S-step
S1 = casos.sossol('S1','mosek',sos1,opts);

% solver 2: V-step
s1sym  = casos.PS.sym('s1',basis(s1));
s2sym  = casos.PS.sym('s2',basis(s2));
s3sym  = casos.PS.sym('s3',basis(s3));
s4sym  = casos.PS.sym('s4',basis(s4));
s51sym = casos.PS.sym('s51',basis(s51));
s52sym = casos.PS.sym('s52',basis(s52));
s61sym = casos.PS.sym('s61',basis(s61));
s62sym = casos.PS.sym('s62',basis(s62));
s8sym = casos.PS.sym('s8',basis(s8));

Vsym = casos.PS.sym('V',basis(V));
ksym = casos.PS.sym('k',basis(k));


sos2 = struct('x',[V;s7;s8],... 
              'p',[ksym;s1sym;s2sym;s3sym;s4sym;s51sym;s52sym;s61sym;s62sym;Vsym;t0;t1;lsym]);   % parameter

sos2.('g') = [s1sym*V   - s2sym*hT - nabla(V,t) - nabla(V,x)*(f + gx*ksym);
              ksym - umin  + s51sym*V - s61sym*hT;
              umax - ksym  + s52sym*V - s62sym*hT;
              s3sym*V   - s4sym*hT - g;
              subs(V,t,t1)    -  lsym]; % ensure we are growing and RS is outside the local terminal set

% states + constraint are SOS cones
opts = struct;
opts.Kx = struct('s', 2,'l',1);
opts.Kc = struct('s', 6);

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

% solver for V-step
S2 = casos.sossol('S2','mosek',sos2,opts);

tend1 = toc;
disp(['took ' num2str(tend1) ' s'])

%% gamma-beta-V-iteration

% initial guess storage function

Vval_old = x'*P*x-1;
l0       = Vval_old;

l = l0;
tic


disp('-----------------------------------------')
disp('Start iteration ... ')
% profile on

for step = 2:N+1

    % set local time interval
    t0 = 0;
    t1 = dt;
        
            % multiplier-step
            sol1 = S1('p',[Vval_old; t0; t1], ...
                      'lbx',klb, ...
                      'ubx',kub);

             switch (S1.stats.UNIFIED_RETURN_STATUS)
                case 'SOLVER_RET_SUCCESS' 
        
                    disp(['S-step feasible in ' num2str(step)])
        
                case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
        
                    disp(['S-step infeasible in ' num2str(step)])
                            
                otherwise, error('Failed.')
            end
        
        
            sol2 = S2('p',[sol1.x; Vval_old; t0; t1;l],'lbx',Vlb,'ubx',Vub);
        
        
            switch (S2.stats.UNIFIED_RETURN_STATUS)
                case 'SOLVER_RET_SUCCESS' 
        
                     Vval = sol2.x(1);
        
                     disp(['V-step feasible in ' num2str(step)])
                    
                case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
        
                    disp(['V-step infeasible in ' num2str(step)])
        
        
                    break
                            
                otherwise, error('Failed.')
            end
         
   
      if step == N+1
            disp('done')
            figure(2)
            pcontour(g,0,[-2 2 -2 2],'k')
            hold on
            pcontour(l0,0,[-2 2 -2 2],'g')
            pcontour((subs(Vval,t,t0)),0,[-2 2 -2 2],'b')
            pcontour(l,0,[-2 2 -2 2],'r')
      else
           figure(2)
            pcontour(g,0,[-2 2 -2 2],'k')
            hold on
            pcontour(l0,0,[-2 2 -2 2],'g')
            pcontour((subs(Vval,t,t0)),0,[-2 2 -2 2],'b')
            pcontour(l,0,[-2 2 -2 2],'r')
           
            % current reachable set becomes new terminal set
            l        = subs(Vval,t,t0); % BRS from current time step, becomes terminal set of next
            Vval_old = l;
      end
end
toc

%%
figure(3)
            pcontour(g,0,[-2 2 -2 2],'k')
            hold on
            pcontour(l0,0,[-2 2 -2 2],'g')
            pcontour((subs(Vval,t,t0)),0,[-2 2 -2 2],'b')
            pcontour(l,0,[-2 2 -2 2],'r')
