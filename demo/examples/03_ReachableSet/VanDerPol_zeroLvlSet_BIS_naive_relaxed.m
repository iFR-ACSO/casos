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
T  = 5;
dt = 1;

N = T/dt;

maxIter = N;

% time polynomial
hT = (t-t0)*(t1-t);

% hT = t/t1;
% Storage function and control law

Vmonom = monomials([x;t],0:4);
% Vmonom(5) = [];
Vmonom(5) = [];
Vmonom(10) = [];
Vmonom(19) = [];


V = casos.PS.sym('v',Vmonom);
k = casos.PS.sym('k',monomials([x;t],0:4));

lsym  = casos.PS.sym('l',basis(subs(V,t,t1)));

% SOS multiplier

s1 = casos.PS.sym('s1',monomials([x;t],0:4),'gram'); % dissipation ineq.
s2 = casos.PS.sym('s2',monomials([x;t],4),'gram');  % dissipation ineq.

s51 = casos.PS.sym('s51',monomials([x;t],0:2),'gram'); % control constraint
s52 = casos.PS.sym('s52',monomials([x;t],0:2),'gram'); % control constraint

s61 = casos.PS.sym('s61',monomials([x;t],2),'gram');   % control constraint
s62 = casos.PS.sym('s62',monomials([x;t],2),'gram');   % control constraint

s7 = casos.PS.sym('s7',monomials(x,0:2),'gram');   % grow constraint
s8 = casos.PS.sym('s8',monomials(x,0:4),'gram');   % grow constraint


%% setup solver

% solver 1: sos-multiplier and control law
tic
disp('-----------------------------------------')
disp('Building solver ... ')

sos1 = struct('x',[k;s1;s2;s51;s52;s61;s62;s8],'p',[V;t0;t1;lsym]);     % parameter

sos1.('g') = [s1*V - s2*hT - nabla(V,t) - nabla(V,x)*(f + gx*k);
              k - umin   + s51*V - s61*hT;
              umax - k   + s52*V - s62*hT;
              s8*subs(V,t,t1) - subs(V,t,t0)];

% states + constraint are SOS cones
opts    = struct;
opts.Kx = struct('s', 7,'l',1);
opts.Kc = struct('s', 4);
opts.error_on_fail = 0;
klb = casos.PS(basis(k),-inf);
kub = casos.PS(basis(k),+inf);

% solver for S-step
S1 = casos.sossol('S1','mosek',sos1,opts);

% solver 2: V-step
s1sym  = casos.PS.sym('s1',basis(s1));
s2sym  = casos.PS.sym('s2',basis(s2));
s51sym = casos.PS.sym('s51',basis(s51));
s52sym = casos.PS.sym('s52',basis(s52));
s61sym = casos.PS.sym('s61',basis(s61));
s62sym = casos.PS.sym('s62',basis(s62));
s7sym = casos.PS.sym('s7',basis(s7));
s8sym = casos.PS.sym('s8',basis(s8));

Vsym = casos.PS.sym('V',basis(V));
ksym = casos.PS.sym('k',basis(k));


sos2 = struct('x',[V;s7],...
              'p',[ksym;s1sym;s2sym;s51sym;s52sym;s61sym;s62sym;s8sym;Vsym;t0;t1;lsym]);   % parameter

sos2.('g') = [s1sym*V   - s2sym*hT - nabla(V,t) - nabla(V,x)*(f + gx*ksym);
              ksym - umin  + s51sym*V - s61sym*hT;
              umax - ksym  + s52sym*V - s62sym*hT;
              subs(V,t,t1)    -  s7*lsym;               % compare to Yin; lsym has the opposite sign
              s8sym*subs(Vsym,t,t1) - subs(V,t,t0)    ]; % ensure we are growing and RS is outside the local terminal set

% states + constraint are SOS cones
opts = struct;
opts.Kx = struct('s', 1,'l',1);
opts.Kc = struct('s', 5);

opts.error_on_fail = 0;
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

N = 100;

% for step = 1:N %2:N+1
iter = 0;
while 1

    % set local time interval
    t0 = 0; %T-(step-1)*dt;
    t1 = 1; %T-(step-2)*dt;
      

    % dt = t1-t0;

            % multiplier-step
            sol1 = S1('p',[Vval_old; t0; t1;l], ...
                      'lbx',klb, ...
                      'ubx',kub);

             switch (S1.stats.UNIFIED_RETURN_STATUS)
                case 'SOLVER_RET_SUCCESS' 
        
                    disp(['S-step feasible in ' num2str(iter)])
        
                case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
        
                    disp(['S-step infeasible in ' num2str(iter)])
                            
                  return
            end
        
        
            sol2 = S2('p',[sol1.x; Vval_old; t0; t1;l],'lbx',Vlb,'ubx',Vub);
        
        
            switch (S2.stats.UNIFIED_RETURN_STATUS)
                case 'SOLVER_RET_SUCCESS' 
        
                     Vval = sol2.x(1);
                     Vval = cleanpoly(Vval,1e-12);
                     disp(['V-step feasible in ' num2str(iter)])
                    
                case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}
        
                    disp(['V-step infeasible in ' num2str(iter)])
        
        
                    return
                            
                % otherwise, error('Failed.')
            end
         
   
      % if iter == N+1
      %       disp('done')
      %       figure(2)
      %       clf
      %       % pcontour(g,0,[-2 2 -2 2],'k')
      %       hold on
      %       pcontour(l0,0,[-2 2 -2 2],'g')
      %       pcontour((subs(Vval,t,t0)),0,[-2 2 -2 2],'b')
      %       % pcontour(l,0,[-2 2 -2 2],'r')
      % else
           figure(2)
           clf
            pcontour(g,0,[-2 2 -2 2],'k')
            hold on
            pcontour(l0,0,[-2 2 -2 2],'g')
            pcontour((subs(Vval,t,t0)),0,[-2 2 -2 2],'b')
            pcontour(l,0,[-2 2 -2 2],'r')
             % Vval
            % current reachable set becomes new terminal set
            l        = subs(Vval,t,t0); % BRS from current time step, becomes terminal set of next
            Vval_old = l;
      % end


      iter  = iter +1;
end
toc

% %%
% figure(3)
%             % pcontour(g,0,[-2 2 -2 2],'k')
%             hold on
%             pcontour(l0,0,[-2 2 -2 2],'g')
%             pcontour((subs(Vval,t,t0)),0,[-2 2 -2 2],'b')
%             pcontour(l,0,[-2 2 -2 2],'r')
