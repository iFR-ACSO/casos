% Inner-approximation reachable set Van-der-Pol Oscillator

clear
close all
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

% time polynomial
hT = (t)*(T-t);

% Storage function and control law
V = casos.PS.sym('v',monomials([x;t],0:2));
k = casos.PS.sym('k',monomials([t;x],0:4));

% SOS multiplier

s1 = casos.PS.sym('s1',monomials([x;t],0:4),'gram'); % dissipation ineq.
s2 = casos.PS.sym('s2',monomials([x;t],4),'gram');  % dissipation ineq.

s3 = casos.PS.sym('s3',monomials([x;t],0:2),'gram');% state constraint
s4 = casos.PS.sym('s4',monomials([x;t],2),'gram');  % state constraint

s51 = casos.PS.sym('s51',monomials([x;t],0:2),'gram'); % control constraint
s52 = casos.PS.sym('s52',monomials([x;t],0:2),'gram'); % control constraint

s61 = casos.PS.sym('s61',monomials([x;t],2),'gram');   % control constraint
s62 = casos.PS.sym('s62',monomials([x;t],2),'gram');   % control constraint

s7 = casos.PS.sym('s7',monomials(x,0:2),'gram');   % grow constraint


cost = 0;
N = 100;
counter = 0;
while true

a = -1; 
b = 1; 
r = (b-a).*rand(2,1) + a; 

        
        if double(subs(g,x,r)) <= 0
            plot(r(1),r(2),'k*')
            hold on
    
            cost = cost + subs(subs(V,t,0),x,r);
    
            counter = counter +1;
            if counter >= N
                break
            end
       end
end

pcontour(g,0,[-3 3 -3 3],'r')

cost = cost*1/N;


%% setup solver

% solver 1: sos-multiplier and control law
tic
disp('-----------------------------------------')
disp('Building solver ... ')

% profile on
sos1 = struct('x',[k,s1,s2,s3,s4,s51,s52,s61,s62]','p',V);     % parameter

sos1.('g') = [s1*V - s2*hT - nabla(V,t) - nabla(V,x)*(f + gx*k);
              k - umin  + s51*V - s61*hT;
              umax - k   + s52*V - s62*hT;
              s3*V - s4*hT - g];

% states + constraint are SOS cones
opts    = struct;
opts.error_on_fail = 0;
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

Vsym = casos.PS.sym('V',basis(V));
ksym = casos.PS.sym('k',basis(k));


sos2 = struct('x',V',...                                                                % dec. variable
              'f',cost,...
              'p',[s1sym,s2sym,s3sym,s4sym,s51sym,s52sym,s61sym,s62sym,ksym,Vsym,t0,t1]');   % parameter

sos2.('g') = [s1sym*V   - s2sym*hT - nabla(V,t) - nabla(V,x)*(f + gx*ksym);
              ksym - umin  + s51sym*V - s61sym*hT;
              umax - ksym  + s52sym*V - s62sym*hT;
              s3sym*V   - s4sym*hT - g;
              subs(V,t,t1)    -  l];

% states + constraint are SOS cones
opts = struct;
opts.error_on_fail = 0;
opts.Kx = struct('s', 0,'l',1);
opts.Kc = struct('s', 5);

Vlb = casos.PS(basis(V),-inf);
Vub = casos.PS(basis(V),+inf);

% solver for V-step
S2 = casos.sossol('S2','mosek',sos2,opts);


tend1 = toc;
disp(['took ' num2str(tend1) ' s'])
% profile viewer
% profile off

%% gamma-beta-V-iteration

% initial guess storage function
Vval = x'*P*x-1;
tic

disp('-----------------------------------------')
disp('Start iteration ... ')
% profile on
for iter = 1:100

    % multiplier-step
    sol1 = S1('p',Vval,'lbx',klb,'ubx',kub);
     switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

            kval   = sol1.x(1);
            s1val  = sol1.x(2);
            s2val  = sol1.x(3);
            s3val  = sol1.x(4);
            s4val  = sol1.x(5);
            s51val = sol1.x(6);
            s52val = sol1.x(7);
            s61val = sol1.x(8);
            s62val = sol1.x(9);

            disp(['S-step feasible in ' num2str(iter)])

        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter)])
                    
         otherwise
             break
    end


    sol2 = S2('p',[s1val,s2val,s3val,s4val,s51val,s52val,s61val,s62val,kval,Vval,0,T],'lbx',Vlb,'ubx',Vub);


    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 

             Vval = sol2.x(1);

             disp(['V-step feasible in ' num2str(iter)])
            
        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['V-step infeasible in ' num2str(iter)])

                break
        otherwise
            disp(['V-step infeasible in ' num2str(iter)])
            break
    end
  

    if iter  < 2
       oldcost = sol2.f
    else
        if abs(double(sol2.f-oldcost)) < 1e-3
            disp('Converged')
            break
        end
        oldcost = sol2.f
    end

end

tendIter = toc;


figure()
pcontour((subs(Vval,t,0)),0,[-1 1 -1 1])
hold on
for k = 0.1:0.1:T-0.1
   pcontour((subs(Vval,t,k)),0,[-1 1 -1 1],'g--') 
end
pcontour((subs(Vval,t,T)),0,[-1 1 -1 1],'b--')
pcontour(g,0,[-1 1 -1 1],'k-')
pcontour(l,0,[-1 1 -1 1],'r-')
% 
%  domain = [x(1), -1, 1; x(2), -1, 1];
% volume = pvolume(subs(Vval,t,0),0,domain)