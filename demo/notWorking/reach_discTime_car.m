clear
close all
clc

%% define variables
x = casos.PS('x',2,1);
u = casos.PS('u',1);

Ts = 20e-3; % sampling time
[Ad,Bd,ulim,xlim] = car2d(Ts);

N = 1;

% scaling
D = 2*diag([1.5 3])^-1;
d = .75;

% open system
f = @(x,u) D\Ad*D*x + D\Bd*d*u;

%% Terminal Penalty & Constraints
% terminal LQR
[k,P] = dlqr(Ad, Bd, 30*eye(2), 1);

% terminal set constraint
P0 = .2;
p0 = x'*D*P*D*x - P0;

umin = d\ulim(1);
umax = d\ulim(2);

xlim = D\xlim/3;

g = (sum((x./xlim(:,2)).^2) - 1);

q  = x'*D*x;
h0 = -k*D*x;

% admissible inputs
hu =  -[u-umin; 
        umax-u];

% initial guess for control laws and storage functions
Vval = p0*ones([N+1 1]);
hval = h0*ones([N 1]);


%% solver 1
t  = casos.PS.sym('t',monomials([x;u],1:2),[1 N]);
s2 = casos.PS.sym('s2',monomials(x,1:3),[2 N],'gram');
s0 = casos.PS.sym('s0',monomials(x,1:2),[2 N],'gram');
sg = casos.PS.sym('sg',monomials(x,0:1),[length(g) N],'gram');

V = casos.PS.sym('V',monomials(x,0:2),[N+1 1]);
h = casos.PS.sym('h',monomials(x,1),[N 1]);

Vsym = casos.PS.sym('V',basis(V(1)),[N+1 1]);
hsym = casos.PS.sym('h',basis(h(1)),[N 1]);

% setup constraint function
constraints = [];

for i = 1:N
    constraints =  [constraints;   
                    -( t(i)*(u - hsym(i)) + subs(Vsym(i+1),x,f(x,u)) - Vsym(i) );
                    -( -s2(:,i)*Vsym(i) + (s0(:,i)+1).*subs(hu,u,hsym(i)) );
                     sg(:,i)*Vsym(i) - g];
end


sos1 = struct('x',[t;s2;s0;sg], ...
              'g', constraints,...
              'p',[Vsym ;hsym ]);     % parameter


opts                = struct;
opts.error_on_fail  = 1;
opts.Kx             = struct('s', length(sos1.x)-length(t),'l',length(t));
opts.Kc = struct('s', length(sos1.g));

% bounds for gram-like variables
tlb = casos.PS(basis(t),-inf);
tub = casos.PS(basis(t),+inf);

% solver for s-step
S1 = casos.sossol('S1','mosek',sos1,opts);


%% solver 2
sb = casos.PS.sym('sb',monomials(x,0:1),'gram');
V0 = casos.PS.sym('V0',monomials(x,0:2));


tsym  = casos.PS.sym('t',basis(t(1)),[1 N]);
s2sym = casos.PS.sym('s2',basis(s2(1)),[2 N]);
s0sym = casos.PS.sym('s0',basis(s0(1)),[2 N]);
sgsym = casos.PS.sym('sg',basis(sg(1)),[length(g) N]);

constraints = [];
for i = 1:N

    constraints = [constraints; 
                    -( tsym(i)*(u - h(i)) + subs(V(i+1),x,f(x,u)) - V(i)  );
                    -( -s2sym(:,i)*V(i) + (s0sym(:,i)+1).*subs(hu,u,h(i)) );
                     sgsym(:,i)*V(i) - g];
end


constraints = [constraints; 
                V(end) - p0; 
                -(-sb*Vsym(1) + V(1) )];


sos2 = struct('x',[V;h;sb], ...
              'g', constraints,...
              'p',[Vsym;tsym;s2sym;s0sym;sgsym]);     % parameter


opts = struct;
opts.error_on_fail = 1;
opts.Kx = struct('s', length(sb), ...
                 'l', length(V) + length(h));
opts.Kc = struct('s', length(sos2.g));


Vlb = casos.PS(basis(V),-inf,[N+1 1]);
Vub = casos.PS(basis(V),+inf,[N+1 1]);

hlb = casos.PS(basis(h),-inf,[N 1]);
hub = casos.PS(basis(h),+inf,[N 1]);

% solver for V-step
S2 = casos.sossol('S2','mosek',sos2,opts);


for iter = 1:1

    sol1 = S1('p',[Vval;hval],'lbx',tlb,'ubx',tub);

    switch (S1.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 
                    
            disp(['S-step feasible in ' num2str(iter)])



        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter)])

        otherwise
    end


    sol2 = S2('p',[Vval;sol1.x], ...
              'lbx',[Vlb;hlb], ...
              'ubx',[Vub;hub]);


    switch (S2.stats.UNIFIED_RETURN_STATUS)
        case 'SOLVER_RET_SUCCESS' 
                    

            disp(['V-step feasible in ' num2str(iter)])
            
            Vval = sol2.x(1:2)
            hval = sol2.x(3)


        case {'SOLVER_RET_INFEASIBLE' 'SOLVER_RET_NAN'}

            disp(['S-step infeasible in ' num2str(iter)])

        otherwise
    end






end